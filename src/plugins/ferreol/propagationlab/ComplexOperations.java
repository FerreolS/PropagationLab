/**
 *
 */
package plugins.ferreol.propagationlab;

import icy.sequence.Sequence;
import mitiv.array.DoubleArray;
import mitiv.array.ShapedArray;
import mitiv.base.Shape;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.mitiv.io.Icy2TiPi;

/**
 * @author ferreol
 *
 */
public class ComplexOperations extends EzPlug implements Block, EzStoppable {

    protected EzVarSequence input1;
    protected EzVarText       inputOption1;
    protected EzVarSequence input2;
    protected EzVarText       inputOption2;

    protected EzVarSequence output;
    protected Sequence outputSequence;


    protected EzVarText       outputOption;  // Combobox for variance estimation
    protected final static String[] outputOptions = new String[]{"Cartesian","Polar","Real part","Imaginary part","modulus","phase"};

    @Override
    protected void initialize()
    {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        input1 = new EzVarSequence("input 1");
        // inputOption1 = new EzVarText(      "type:", outputOptions, false);

        input2 = new EzVarSequence("input 2");
        inputOption2 = new EzVarText(      "type:", outputOptions, false);

        outputOption = new EzVarText(      "Output:", outputOptions, false);

        if (isHeadLess()) {
            output = new EzVarSequence("Propagated field Image");
        }
    }
    @Override
    protected void execute() {
        Sequence inputSequence1 = input1.getValue();
        ShapedArray inputArray1 = Icy2TiPi.sequenceToArray(inputSequence1);
        Sequence inputSequence2 = input2.getValue();
        ShapedArray inputArray2 = Icy2TiPi.sequenceToArray(inputSequence1);
        Sequence outputSequence= new Sequence();
        outputSequence.copyMetaDataFrom(inputSequence1, false);
        DoubleArray outputArray;

        Shape outputShape = inputArray1.getShape();

        if (!outputShape.equals(inputArray2.getShape())){
            throw new IllegalArgumentException("Both input should have the same shape");
        }

        outputArray.as1D().toDouble().flatten()

    }
}