/**
 *
 */
package plugins.ferreol.PropagationLab;

import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import icy.util.StringUtil;
import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.Double4D;
import mitiv.array.ShapedArray;
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
public class ComplexConvert extends EzPlug implements Block, EzStoppable,PluginBundled {
    protected EzVarSequence input;
    protected Sequence inputSequence;
    //  protected ShapedArray inputArray;
    protected ShapedArray outputArray;
    protected EzVarSequence output;
    protected Sequence outputSequence;

    protected EzVarText       inputOption;  // Combobox
    protected final static String[] inputOptions = new String[]{"Cartesian","Polar"};


    protected EzVarText       outputOption;  // Combobox
    protected final static String[] outputOptions = new String[]{"Cartesian","Polar","Real part","Imaginary part","modulus","phase","Log(modulus)"};

    @Override
    protected void initialize()
    {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        input = new EzVarSequence("input");
        outputOption = new EzVarText(      "Output type:", outputOptions, false);
        inputOption = new EzVarText(      "input  type:", inputOptions, false);

        addEzComponent(input);
        addEzComponent(inputOption);
        addEzComponent(outputOption);

        if (isHeadLess()) {
            output = new EzVarSequence("Output");
        }
    }
    @Override
    protected void execute() {
        Sequence inputSequence = input.getValue();
        if (inputSequence==null){
            return;
        }
        outputArray = Icy2TiPi.sequenceToArray(inputSequence);
        Sequence outputSequence= new Sequence();
        outputSequence.copyMetaDataFrom(inputSequence, false);


        if(StringUtil.equals(inputOption.getValue(),outputOption.getValue() )){


            if (isHeadLess()) {
                output.setValue(outputSequence);
            }
            return;
        }
        int complexdim = Math.min(outputArray.getRank()-1,2);
        int rank = outputArray.getRank();

        if((rank<4)&&(outputArray.getDimension(complexdim)==2)){ //Assuming complex input
            outputArray = outputArray.movedims( complexdim,0).toDouble().copy();
        }else{
            throw new IllegalArgumentException("Input must be complex cartesian with less than 4D");
        }


        if(StringUtil.equals(inputOption.getValue(),inputOptions[0])){ // Cartesian input
            if(StringUtil.equals(outputOption.getValue(),outputOptions[2]) ){//Real part
                switch (outputArray.getRank()){
                    case 2:
                        IcyImager.show( ((Double2D) outputArray).slice(0,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 3:
                        IcyImager.show( ((Double3D) outputArray).slice(0,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 4:
                        IcyImager.show( ((Double4D) outputArray).slice(0,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                }
                outputSequence.setChannelName(0, "Real part");
            }else if(StringUtil.equals(outputOption.getValue(),outputOptions[3]) ){// imaginary part
                switch (outputArray.getRank()){
                    case 2:
                        IcyImager.show( ((Double2D) outputArray).slice(1,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 3:
                        IcyImager.show( ((Double3D) outputArray).slice(1,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 4:
                        IcyImager.show( ((Double4D) outputArray).slice(1,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                }
                outputSequence.setChannelName(0, "Imaginary part");
            }else{
                if(StringUtil.equals(outputOption.getValue(),outputOptions[1] )){ ///Polar output
                    for(int i=0;i<outputArray.getNumber();i=i+2){
                        double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                        double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                        ((Double1D)outputArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
                        ((Double1D)outputArray.toDouble().as1D()).set(i+1,Math.atan2(im,re));
                    }
                    IcyImager.show(outputArray,outputSequence,0,inputSequence.getName(), isHeadLess() );
                } else{
                    if(StringUtil.equals(outputOption.getValue(),outputOptions[4] )){//modulus

                        for(int i=0;i<outputArray.getNumber();i=i+2){
                            double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                            double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                            ((Double1D)outputArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
                        }

                    }else if(StringUtil.equals(outputOption.getValue(),outputOptions[5])){//phase

                        for(int i=0;i<outputArray.getNumber();i=i+2){
                            double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                            double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                            ((Double1D)outputArray.toDouble().as1D()).set(i,Math.atan2(im,re));
                        }

                    }else if(StringUtil.equals(outputOption.getValue(),outputOptions[6] )){//log modulus

                        for(int i=0;i<outputArray.getNumber();i=i+2){
                            double re = ((Double1D)outputArray.toDouble().as1D()).get(i);
                            double im = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                            ((Double1D)outputArray.toDouble().as1D()).set(i, Math.log10(re*re+im*im+1E-15));
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
                }
                if(StringUtil.equals(outputOption.getValue(),outputOptions[5] )){//phase
                    outputSequence.setChannelName(0, "Phase");

                }else if(StringUtil.equals(outputOption.getValue(),outputOptions[6]) ){//squared modulus
                    outputSequence.setChannelName(0, "Log(modulus)");
                }else{
                    outputSequence.setChannelName(0, "Modulus");
                    if(StringUtil.equals(outputOption.getValue(),outputOptions[1])){
                        outputSequence.setChannelName(1, "Phase");
                    }
                }

            }
        }else{ // Polar input

            for(int i=0;i<outputArray.getNumber();i=i+2){
                double mod = ((Double1D)outputArray.toDouble().as1D()).get(i);
                double ph = ((Double1D)outputArray.toDouble().as1D()).get(i+1);
                ((Double1D)outputArray.toDouble().as1D()).set(i,mod*Math.cos(ph));
                ((Double1D)outputArray.toDouble().as1D()).set(i+1,mod*Math.sin(ph));
            }

            if(StringUtil.equals(outputOption.getValue(),outputOptions[0] )){
                IcyImager.show(outputArray,outputSequence,0,inputSequence.getName(), isHeadLess() );
                outputSequence.setChannelName(0, "Real part");
                outputSequence.setChannelName(1, "Imaginary part");
            }
            if(StringUtil.equals(outputOption.getValue(),outputOptions[2]) ){//Real part
                switch (outputArray.getRank()){
                    case 2:
                        IcyImager.show( ((Double2D) outputArray).slice(0,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 3:
                        IcyImager.show( ((Double3D) outputArray).slice(0,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 4:
                        IcyImager.show( ((Double4D) outputArray).slice(0,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                }
                outputSequence.setChannelName(0, "Real part");
            }else if(StringUtil.equals(outputOption.getValue(),outputOptions[3]) ){// imaginary part
                switch (outputArray.getRank()){
                    case 2:
                        IcyImager.show( ((Double2D) outputArray).slice(1,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 3:
                        IcyImager.show( ((Double3D) outputArray).slice(1,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                    case 4:
                        IcyImager.show( ((Double4D) outputArray).slice(1,0),outputSequence,inputSequence.getName(), isHeadLess() );
                        break;
                }
                outputSequence.setChannelName(0, "Imaginary part");
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
        inputMap.add("input type", inputOption.getVariable());
        inputMap.add("output type", outputOption.getVariable());
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
    /* (non-Javadoc)
     * @see icy.plugin.interface_.PluginBundled#getMainPluginClassName()
     */
    @Override
    public String getMainPluginClassName() {
        return "PropagationLab";
    }
}
