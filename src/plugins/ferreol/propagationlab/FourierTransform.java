
package plugins.ferreol.propagationlab;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_3D;
import icy.sequence.Sequence;
import mitiv.array.ArrayFactory;
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

public class FourierTransform  extends EzPlug  implements Block, EzStoppable{
    protected EzVarSequence input;
    protected Sequence inputSequence;
    protected ShapedArray inputArray;
    protected DoubleArray outputArray;
    protected EzVarSequence output;
    protected Sequence outputSequence;
    protected EzVarBoolean    direction;

    protected EzVarText       outputOption;  // Combobox for variance estimation
    protected final String[] outpoutOptions = new String[]{"Cartesian","Polar","Real part","Imaginary part","modulus","phase","Squared modulus"};
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
        outputOption = new EzVarText(      "Output:", outpoutOptions, false);
        addEzComponent(input);
        addEzComponent(direction);
        addEzComponent(outputOption);
        if (isHeadLess()) {
            output = new EzVarSequence("Output Image");
        }
    }
    @Override
    protected void execute() {
        Sequence inputSequence = input.getValue();
        inputArray = Icy2TiPi.sequenceToArray(inputSequence);

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


        IcyImager.show(outputArray,outputSequence,0,"Fourier transform of "+inputSequence.getName(), isHeadLess() );

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
        // TODO Auto-generated method stub

    }
    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        // TODO Auto-generated method stub

    }

}
