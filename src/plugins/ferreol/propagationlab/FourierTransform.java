
package plugins.ferreol.propagationlab;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.ArrayFactory;
import mitiv.array.ArrayUtils;
import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.Double4D;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.Icy2TiPi;
import plugins.mitiv.io.IcyImager;

public class FourierTransform  extends EzPlug  implements Block, EzStoppable, PluginBundled{
    protected EzVarSequence input;
    protected Sequence inputSequence;
    protected ShapedArray inputArray;
    protected DoubleArray outputArray;
    protected EzVarSequence output;
    protected Sequence outputSequence;
    protected EzVarBoolean    direction;
    protected EzVarBoolean    fftshift;


    protected EzVarText       outputOption;  // Combobox for variance estimation
    protected final static String[] outputOptions = new String[]{"Cartesian","Polar","Real part","Imaginary part","modulus","phase","Squared modulus"};
    private boolean iscomplex;
    private Shape outputShape;



    @Override
    protected void initialize()
    {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        input = new EzVarSequence("input");
        input.setToolTipText("input with real and imaginary part in channel 0 and 1 respectively");
        direction = new EzVarBoolean("Backward", false);
        direction.setToolTipText("Direction of the transform (backward if checked");
        fftshift = new EzVarBoolean("FFT shift", false);
        fftshift.setToolTipText("Swap quadrant to center the 0 frequency");
        outputOption = new EzVarText(      "Output:", outputOptions, false);
        addEzComponent(input);
        addEzComponent(direction);
        addEzComponent(fftshift);
        addEzComponent(outputOption);
        if (isHeadLess()) {
            output = new EzVarSequence("Output Image");
        }
    }
    @Override
    protected void execute() {
        Sequence inputSequence = input.getValue();
        inputArray = Icy2TiPi.sequenceToArray(inputSequence);
        Sequence outputSequence= new Sequence();
        outputSequence.copyMetaDataFrom(inputSequence, false);
        double[] data;

        if(((inputArray.getRank()<3)&&(inputArray.getDimension(inputArray.getRank()-1)==2)) ||((inputArray.getRank()>2) &&(inputArray.getDimension(2)==2))){ //Assuming complex input
            iscomplex = true;
            inputArray = inputArray.movedims(2, 0);
            data  = inputArray.toDouble().flatten();
            switch (inputArray.getRank()) {
                case 2:{
                    DoubleFFT_1D FFT = new DoubleFFT_1D(inputArray.getDimension(1));

                    if (direction.getValue()){
                        FFT.complexInverse(data, true);
                    }else{
                        FFT.complexForward(data);
                    }
                }
                break;

                case 3:{
                    DoubleFFT_2D FFT = new DoubleFFT_2D(inputArray.getDimension(1), inputArray.getDimension(2));

                    if (direction.getValue()){
                        FFT.complexInverse(data, true);
                    }else{
                        FFT.complexForward(data);
                    }
                }
                break;

                case 4:{
                    DoubleFFT_3D FFT = new DoubleFFT_3D(inputArray.getDimension(1), inputArray.getDimension(2), inputArray.getDimension(3));

                    if (direction.getValue()){
                        FFT.complexInverse(data, true);
                    }else{
                        FFT.complexForward(data);
                    }
                    break;
                }
                default:
                    break;
            }
            outputShape=inputArray.getShape();
        }else{ // Assuming a real input
            iscomplex = false;
            int [] newdims = new  int[inputArray.getRank()+1];
            newdims[0] = 2;
            for (int i = 1; i < newdims.length; i++) {
                newdims[i] = inputArray.getDimension(i-1);
            }
            data = java.util.Arrays.copyOf(inputArray.toDouble().flatten(),2*inputArray.getNumber());

            switch (inputArray.getRank()) {
                case 1:{
                    DoubleFFT_1D FFT = new DoubleFFT_1D(inputArray.getDimension(0));

                    if (direction.getValue()){
                        FFT.realInverseFull(data, true);
                    }else{
                        FFT.realForwardFull(data);
                    }
                }
                break;

                case 2:{
                    DoubleFFT_2D FFT = new DoubleFFT_2D(inputArray.getDimension(0), inputArray.getDimension(1));

                    if (direction.getValue()){
                        FFT.realInverseFull(data, true);
                    }else{
                        FFT.realForwardFull(data);
                    }
                }
                break;

                case 3:{
                    DoubleFFT_3D FFT = new DoubleFFT_3D(inputArray.getDimension(0), inputArray.getDimension(1), inputArray.getDimension(2));

                    if (direction.getValue()){
                        FFT.realInverseFull(data, true);
                    }else{
                        FFT.realForwardFull(data);
                    }
                }
                default:
                    break;
            }
            outputShape = new Shape(newdims);
        }
        outputArray = ArrayFactory.wrap(data, outputShape);

        if (fftshift.getValue()){
            int rank = outputArray.getRank();
            int[] off = new int[rank];
            for (int k = 1; k < rank; ++k) {
                int dim = outputArray.getDimension(k);
                off[k] = -(dim/2);
            }

            outputArray = (DoubleArray) ArrayUtils.roll(outputArray,off).copy();

        }

        if(outputOption.getValue()==outputOptions[0] ){ // Cartesian
            IcyImager.show(outputArray,outputSequence,0,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
            outputSequence.setChannelName(0, "Real part");
            outputSequence.setChannelName(1, "Imaginary part");

        }else if(outputOption.getValue()==outputOptions[2] ){//Real part
            switch (outputArray.getRank()){
                case 2:
                    IcyImager.show( ((Double2D) outputArray).slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                    break;
                case 3:
                    IcyImager.show( ((Double3D) outputArray).slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                    break;
                case 4:
                    IcyImager.show( ((Double4D) outputArray).slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                    break;
            }
            outputSequence.setChannelName(0, "Real part");
        }else if(outputOption.getValue()==outputOptions[3] ){// imaginary part
            switch (outputArray.getRank()){
                case 2:
                    IcyImager.show( ((Double2D) outputArray).slice(1,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                    break;
                case 3:
                    IcyImager.show( ((Double3D) outputArray).slice(1,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                    break;
                case 4:
                    IcyImager.show( ((Double4D) outputArray).slice(1,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                    break;
            }
            outputSequence.setChannelName(0, "Imaginary part");
        }else{
            if(outputOption.getValue()==outputOptions[1] ){ //Polar

                for(int i=0;i<outputArray.getNumber();i=i+2){
                    double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                    double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                    ((Double1D)outputArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
                    ((Double1D)outputArray.toDouble().as1D()).set(i+1,Math.atan2(im,re));
                }
                IcyImager.show(outputArray,outputSequence,0,"Fourier transform of "+inputSequence.getName(), isHeadLess() );

            }else{
                if(outputOption.getValue()==outputOptions[4] ){//modulus

                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
                    }

                }else if(outputOption.getValue()==outputOptions[5] ){//phase

                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,Math.atan2(im,re));
                    }

                }else if(outputOption.getValue()==outputOptions[6] ){//squared modulus

                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,(re*re+im*im));
                    }
                }
                switch (outputArray.getRank()){
                    case 2:
                        IcyImager.show( ((Double2D) outputArray).slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                        break;
                    case 3:
                        IcyImager.show( ((Double3D) outputArray).slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                        break;
                    case 4:
                        IcyImager.show( ((Double4D) outputArray).slice(0,0),outputSequence,"Fourier transform of "+inputSequence.getName(), isHeadLess() );
                        break;

                }
                //  outputSequence.getFirstViewer().getLut().getLutChannel(0).setColorMap(new IceColorMap(),false);
            }
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


        // outputSequence.getFirstViewer().getLut().getLutChannel(0).setColorMap(new IceColorMap(),false);


        if (isHeadLess()) {
            output.setValue(outputSequence);
        }
    }

    @Override
    public void clean() {
        // TODO Auto-generated by Icy4Eclipse
    }
    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareInput(VarList inputMap) {

        initialize();
        inputMap.add("input", input.getVariable());
        inputMap.add("direction", direction.getVariable());
        inputMap.add("fftshift", fftshift.getVariable());
        inputMap.add("outputOption", outputOption.getVariable());

    }
    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("output", output.getVariable());
    }
    /* (non-Javadoc)
     * @see icy.plugin.interface_.PluginBundled#getMainPluginClassName()
     */
    @Override
    public String getMainPluginClassName() {
        // TODO Auto-generated method stub
        return null;
    }

}
