/**
 *
 */
package plugins.ferreol.PropagationLab;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.Double3D;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import mitiv.utils.MathUtils;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.ferreol.demics.ToolTipText;
import plugins.mitiv.io.Icy2TiPi;
import plugins.mitiv.io.IcyImager;

/**
 * @author ferreol
 *
 */
public class Propagate extends EzPlug implements Block, EzStoppable, PluginBundled {


    // optical parameters
    protected EzVarDouble     dxy_nm;  //  pixels size in (x,y) in nm
    protected EzVarDouble     lambda;         //  wavelength
    protected EzVarDouble     ni;             //  refractive index of the immersion index
    protected EzVarDouble     depth;

    protected EzVarSequence input;
    protected Sequence inputSequence;
    protected ShapedArray inputArray;
    protected Double3D outputArray;
    protected EzVarSequence output;
    protected Sequence outputSequence;


    protected EzVarText       outputOption;  // Combobox for variance estimation
    protected final static String[] outputOptions = new String[]{"Cartesian","Polar","Real part","Imaginary part","modulus","phase","Squared modulus"};

    int Nxy;

    @Override
    protected void initialize()
    {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        input = new EzVarSequence("input");
        input.setToolTipText("input with real and imaginary part in channel 0 and 1 respectively");
        dxy_nm = new EzVarDouble("dxy(nm):",64.5,0., Double.MAX_VALUE,1.);
        lambda = new EzVarDouble( "\u03BB(nm):",540.,10.,15000.,5);
        ni = new EzVarDouble("refractive index:",1.,0.99,2.,0.05);
        depth = new EzVarDouble("depth (mm):",0,-100,100,.001);
        outputOption = new EzVarText(      "Output:", outputOptions, false);


        addEzComponent(input);
        addEzComponent(dxy_nm);
        dxy_nm.setToolTipText(ToolTipText.doubleDxy);
        addEzComponent(lambda);
        lambda.setToolTipText(ToolTipText.doubleLambda);
        addEzComponent(ni);
        ni.setToolTipText(ToolTipText.doubleNi);
        addEzComponent(depth);
        depth.setToolTipText("Propagation depth");
        addEzComponent(outputOption);

        if (isHeadLess()) {
            output = new EzVarSequence("Propagated field Image");
        }
    }
    @Override
    protected void execute() {
        Sequence inputSequence = input.getValue();
        inputArray = Icy2TiPi.sequenceToArray(inputSequence);
        Sequence outputSequence= new Sequence();
        outputSequence.copyMetaDataFrom(inputSequence, false);

        Nxy = inputArray.getDimension(0);
        if( Nxy != inputArray.getDimension(1)){
            throw new IllegalArgumentException("Can only deal with squared images");
        }

        double pixsize = dxy_nm.getValue()*1E-9 ; //en metres
        double z= depth.getValue()*1E-3;        // en metres
        double wavelenght = lambda.getValue() *1E-9;// en metres
        double k = 2.*Math.PI/wavelenght*ni.getValue();

        double redfreqsize = wavelenght/(pixsize*Nxy*ni.getValue());
        double[] data ;

        DoubleFFT_2D FFT = new DoubleFFT_2D(Nxy, Nxy);
        Shape outputShape;
        if((inputArray.getRank()==3)&&(inputArray.getDimension(2)==2)){ //Assuming complex input
            inputArray = inputArray.movedims(2, 0);
            data  = inputArray.toDouble().flatten();
            FFT.complexForward(data);
            outputShape = inputArray.getShape();
        }else if(inputArray.getRank()==2){ // Assuming real input
            int [] newdims = new  int[3];
            newdims[0] = 2;
            newdims[1] = Nxy;
            newdims[2] = Nxy;
            outputShape = new Shape(newdims);
            data = java.util.Arrays.copyOf(inputArray.toDouble().flatten(),2*inputArray.getNumber());
            FFT.realForwardFull(data);
        }else{
            throw new IllegalArgumentException("Only 2D waves can be propagated");
        }

        //  outputArray =  Double3D.wrap(data,outputShape);

        double r[] = MathUtils.fftDist1D(Nxy, Nxy);

        for (int nx = 0; nx < r.length; nx++) {
            double tmp = Math.pow(r[nx]*redfreqsize,2);
            if (tmp<1){
                double phaseAS = -k*z*Math.sqrt(1 - tmp);
                double cosphase = Math.cos(phaseAS);
                double sinphase = Math.sin(phaseAS);
                double ref= data[2*nx];
                double imf = data[2*nx+1];

                data[2*nx]= ref*cosphase-imf*sinphase;
                data[2*nx+1]=  ref*sinphase+imf*cosphase;
            }else{
                double ev = Math.exp(-k*Math.abs(z)*Math.sqrt(tmp-1));
                data[2*nx]*= ev ;
                data[2*nx+1] *= ev;

            }
        }
        FFT.complexInverse(data, true);
        outputArray =  Double3D.wrap(data,outputShape);

        if(outputOption.getValue()==outputOptions[0] ){ // Cartesian
            IcyImager.show(outputArray,outputSequence,0,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
            outputSequence.setChannelName(0, "Real part");
            outputSequence.setChannelName(1, "Imaginary part");

        }else if(outputOption.getValue()==outputOptions[2] ){//Real part
            IcyImager.show( outputArray.slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
            outputSequence.setChannelName(0, "Real part");
        }else if(outputOption.getValue()==outputOptions[3] ){// imaginary part
            IcyImager.show( outputArray.slice(1,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
            outputSequence.setChannelName(0, "Imaginary part");
        }else{
            if(outputOption.getValue()==outputOptions[1] ){ //Polar
                for(int i=0;i<outputArray.getNumber();i=i+2){
                    double re = outputArray.toDouble().as1D().get(i);
                    double im = outputArray.toDouble().as1D().get(i+1);
                    outputArray.toDouble().as1D().set(i,Math.sqrt(re*re+im*im));
                    outputArray.toDouble().as1D().set(i+1,Math.atan2(im,re));
                }
                IcyImager.show(outputArray,outputSequence,0,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
            }else if(outputOption.getValue()==outputOptions[4] ){//modulus
                for(int i=0;i<outputArray.getNumber();i=i+2){
                    double re = outputArray.toDouble().as1D().get(i);
                    double im = outputArray.toDouble().as1D().get(i+1);
                    outputArray.toDouble().as1D().set(i,Math.sqrt(re*re+im*im));
                }

            }else if(outputOption.getValue()==outputOptions[5] ){//phase

                for(int i=0;i<outputArray.getNumber();i=i+2){
                    double re = outputArray.toDouble().as1D().get(i);
                    double im = outputArray.toDouble().as1D().get(i+1);
                    outputArray.toDouble().as1D().set(i,Math.atan2(im,re));
                }

            }else if(outputOption.getValue()==outputOptions[6] ){//squared modulus

                for(int i=0;i<outputArray.getNumber();i=i+2){
                    double re = outputArray.toDouble().as1D().get(i);
                    double im = outputArray.toDouble().as1D().get(i+1);
                    outputArray.toDouble().as1D().set(i,(re*re+im*im));
                }
            }
            IcyImager.show( outputArray.slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
            if(outputOption.getValue()==outputOptions[5] ){//phase
                outputSequence.setChannelName(0, "Phase");

            }else if(outputOption.getValue()==outputOptions[6] ){//squared modulus
                outputSequence.setChannelName(0, "Squared modulus");
            }else{
                outputSequence.setChannelName(0, "Modulus");
                if(outputOption.getValue()==outputOptions[1]){
                    outputSequence.setChannelName(1, "Phase");
                }
            }
        }


        if (isHeadLess()) {
            output.setValue(outputSequence);
        }
    }


    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareInput(VarList inputMap) {
        initialize();
        inputMap.add("input", input.getVariable());
        inputMap.add("pixel size", dxy_nm.getVariable());
        inputMap.add("wavelength", lambda.getVariable());
        inputMap.add("refractive index", ni.getVariable());
        inputMap.add("depth", depth.getVariable());
        inputMap.add("outputOption", outputOption.getVariable());

    }

    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        // TODO Auto-generated method stub
        outputMap.add("output", output.getVariable());

    }
    /* (non-Javadoc)
     * @see plugins.adufour.ezplug.EzPlug#clean()
     */
    @Override
    public void clean() {
        // TODO Auto-generated method stub

    }

    @Override
    public String getMainPluginClassName() {
        return "PropagationLab";
    }


}
