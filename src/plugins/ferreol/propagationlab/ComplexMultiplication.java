/**
 *
 */
package plugins.ferreol.PropagationLab;

import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import icy.util.StringUtil;
import mitiv.array.ArrayFactory;
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
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.Icy2TiPi;
import plugins.mitiv.io.IcyImager;

/**
 * @author ferreol
 *
 */
public class ComplexMultiplication extends EzPlug implements Block, EzStoppable, PluginBundled {

    protected EzVarSequence input1;
    protected EzVarText       inputOption1;
    protected EzVarSequence input2;
    protected EzVarText       inputOption2;

    protected EzVarSequence output;
    protected Sequence outputSequence;

    Sequence inputSequence1;
    Sequence inputSequence2;


    protected EzVarText       outputOption;  // Combobox for variance estimation
    protected final static String[] outputOptions = new String[]{"Cartesian","Polar","Real part","Imaginary part","modulus","phase","log(modulus)"};

    @Override
    protected void initialize()
    {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        input1 = new EzVarSequence("input 1 (Cartesian)");
        // inputOption1 = new EzVarText(      "type:", outputOptions, false);

        input2 = new EzVarSequence("input 2");
        inputOption2 = new EzVarText(      "input 2 type:", outputOptions, false);

        outputOption = new EzVarText(      "Output type:", outputOptions, false);


        addEzComponent(input1);
        addEzComponent(input2);
        addEzComponent(inputOption2);
        addEzComponent(outputOption);

        if (isHeadLess()) {
            output = new EzVarSequence("Output");
        }
    }
    @Override
    protected void execute() {
        inputSequence1 = input1.getValue();
        inputSequence2 = input2.getValue();


        if ((inputSequence1==null)||(inputSequence2==null)){
            return;
        }
        ShapedArray inputArray1 =  Icy2TiPi.sequenceToArray(inputSequence1);
        ShapedArray inputArray2 = Icy2TiPi.sequenceToArray(inputSequence2);
        Sequence outputSequence= new Sequence();
        //  outputSequence.copyMetaDataFrom(inputSequence1, false);
        DoubleArray outputArray;

        Shape outputShape ;
        int rank = inputArray1.getRank();

        if((rank<4)&&(inputArray1.getDimension(inputArray1.getRank()-1)==2)){ //Assuming complex input
            inputArray1 = inputArray1.movedims( Math.min(inputArray1.getRank()-1,2), 0).toDouble();
            outputShape = inputArray1.getShape();
        }else{
            throw new IllegalArgumentException("First input must be complex cartesian with less than 4D");
        }



        double [] data1 = inputArray1.toDouble().flatten();
        if(StringUtil.equals(inputOption2.getValue(),outputOptions[0]) ){ //Cartesian

            /*   if (!outputShape.equals(inputArray2.getShape())){
                throw new IllegalArgumentException("Both input should have the same shape");
            }*/
            double [] data2 = inputArray2.movedims(2, 0).toDouble().flatten();

            for (int nx = 0; nx < inputArray1.getNumber()/2; nx++){
                double re1= data1[2*nx];
                double im1 = data1[2*nx+1];
                double re2= data2[2*nx];
                double im2 = data2[2*nx+1];
                data1[2*nx]= re1*re2-im1*im2;
                data1[2*nx+1]=  re1*im2+im1*re2;
            }
        }else if(StringUtil.equals(inputOption2.getValue(),outputOptions[1] )){//Polar
            /*        if (!outputShape.equals(inputArray2.getShape())){
                throw new IllegalArgumentException("Both input should have the same shape");
            }
             */
            double [] data2 = inputArray2.movedims(2, 0).toDouble().flatten();

            for (int nx = 0; nx < inputArray1.getNumber()/2; nx++){
                double re1= data1[2*nx];
                double im1 = data1[2*nx+1];
                double mod2= data2[2*nx];
                double phase2 = data2[2*nx+1];
                double cosphase2 = Math.cos(phase2);
                double sinphase2 = Math.sin(phase2);
                data1[2*nx]= mod2*(re1*cosphase2- im1*sinphase2);
                data1[2*nx+1]= mod2*(re1*sinphase2- im1*cosphase2);
            }
        }else if((StringUtil.equals(inputOption2.getValue(),outputOptions[2] ))||(StringUtil.equals(inputOption2.getValue(),outputOptions[4]) )){//Real or modulus
            double [] data2 = inputArray2.toDouble().flatten();

            for (int nx = 0; nx < inputArray1.getNumber()/2; nx++){
                data1[2*nx] *=  data2[nx];
                data1[2*nx+1]*=  data2[nx];
            }

        }else if(StringUtil.equals(inputOption2.getValue(),outputOptions[3] )){//Imaginary
            double [] data2 = inputArray2.toDouble().flatten();

            for (int nx = 0; nx < inputArray1.getNumber()/2; nx++){
                double re1= data1[2*nx];
                double im1 = data1[2*nx+1];
                double im2 = data2[nx];
                data1[2*nx]= -im1*im2;
                data1[2*nx+1]=  re1*im2;
            }

        }else  if(StringUtil.equals(inputOption2.getValue(),outputOptions[5] )){//phase
            double [] data2 = inputArray2.toDouble().flatten();

            for (int nx = 0; nx < inputArray1.getNumber()/2; nx++){
                double re1= data1[2*nx];
                double im1 = data1[2*nx+1];
                double cosphase2 = Math.cos(data2[nx]);
                double sinphase2 = Math.sin(data2[nx]);
                data1[2*nx]= (re1*cosphase2- im1*sinphase2);
                data1[2*nx+1]= (re1*sinphase2- im1*cosphase2);
            }
        }
        outputArray = ArrayFactory.wrap(data1, outputShape);


        if(StringUtil.equals(outputOption.getValue(),outputOptions[0] )){ // Cartesian
            IcyImager.show(outputArray,outputSequence,0,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
            outputSequence.setChannelName(0, "Real part");
            outputSequence.setChannelName(1, "Imaginary part");

        }else if(StringUtil.equals(outputOption.getValue(),outputOptions[2]) ){//Real part
            switch (outputArray.getRank()){
                case 2:
                    IcyImager.show( ((Double2D) outputArray).slice(0,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                    break;
                case 3:
                    IcyImager.show( ((Double3D) outputArray).slice(0,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                    break;
                case 4:
                    IcyImager.show( ((Double4D) outputArray).slice(0,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                    break;
            }
            outputSequence.setChannelName(0, "Real part");
        }else if(StringUtil.equals(outputOption.getValue(),outputOptions[3] )){// imaginary part
            switch (outputArray.getRank()){
                case 2:
                    IcyImager.show( ((Double2D) outputArray).slice(1,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                    break;
                case 3:
                    IcyImager.show( ((Double3D) outputArray).slice(1,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                    break;
                case 4:
                    IcyImager.show( ((Double4D) outputArray).slice(1,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                    break;
            }
            outputSequence.setChannelName(0, "Imaginary part");
        }else{
            if(StringUtil.equals(outputOption.getValue(),outputOptions[1]) ){ //Polar

                for(int i=0;i<outputArray.getNumber();i=i+2){
                    double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                    double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                    ((Double1D)outputArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
                    ((Double1D)outputArray.toDouble().as1D()).set(i+1,Math.atan2(im,re));
                }
                IcyImager.show(outputArray,outputSequence,0,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );

            }else{
                if(StringUtil.equals(outputOption.getValue(),outputOptions[4]) ){//modulus

                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
                    }

                }else     if(StringUtil.equals(outputOption.getValue(),outputOptions[6] )){//modulus

                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,Math.log10(re*re+im*im+1E-15));
                    }

                }else if(StringUtil.equals(outputOption.getValue(),outputOptions[5]) ){//phase

                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,Math.atan2(im,re));
                    }

                }
                switch (outputArray.getRank()){
                    case 2:
                        IcyImager.show( ((Double2D) outputArray).slice(0,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                        break;
                    case 3:
                        IcyImager.show( ((Double3D) outputArray).slice(0,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                        break;
                    case 4:
                        IcyImager.show( ((Double4D) outputArray).slice(0,0),outputSequence,inputSequence1.getName()+ " X " + inputSequence2.getName(), isHeadLess() );
                        break;

                }
                //  outputSequence.getFirstViewer().getLut().getLutChannel(0).setColorMap(new IceColorMap(),false);
            }
            if(StringUtil.equals(outputOption.getValue(),outputOptions[5] )){//phase
                outputSequence.setChannelName(0, "Phase");
            }else{
                outputSequence.setChannelName(0, "Modulus");
                if(StringUtil.equals(outputOption.getValue(),outputOptions[1])){
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
        inputMap.add("input 1", input1.getVariable());
        inputMap.add("input 2", input2.getVariable());
        inputMap.add("input2 type", inputOption2.getVariable());
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
