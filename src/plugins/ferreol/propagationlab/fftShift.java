/**
 *
 */
package plugins.ferreol.PropagationLab;

import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.ArrayUtils;
import mitiv.array.ShapedArray;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.mitiv.io.Icy2TiPi;
import plugins.mitiv.io.IcyImager;

/**
 * @author ferreol
 *
 */
public class fftShift extends EzPlug implements Block, EzStoppable, PluginBundled{
    protected EzVarSequence input;
    protected Sequence inputSequence;
    protected ShapedArray inputArray;
    protected ShapedArray outputArray;
    protected EzVarSequence output;
    protected Sequence outputSequence;

    @Override
    protected void initialize()
    {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        input = new EzVarSequence("input");
        addEzComponent(input);
        if (isHeadLess()) {
            output = new EzVarSequence("Output");
        }
    }
    @Override
    protected void execute() {
        Sequence inputSequence = input.getValue();
        inputArray = Icy2TiPi.sequenceToArray(inputSequence);
        Sequence outputSequence= new Sequence();
        //    outputSequence.copyFrom(inputSequence, false);

        outputSequence.copyMetaDataFrom(inputSequence, false);

        if(((inputArray.getRank()<3)&&(inputArray.getDimension(inputArray.getRank()-1)==2)) ||((inputArray.getRank()>2) &&(inputArray.getDimension(2)==2))){ //Assuming complex input

            int rank = inputArray.getRank();
            int[] off = new int[rank];
            for (int k = 0; k < rank; ++k) {
                int dim = inputArray.getDimension(k);
                off[k] = -(dim/2);
            }
            off[Math.min(inputArray.getRank()-1,2)] = 0;
            inputArray =  ArrayUtils.roll(inputArray,off).copy();
            IcyImager.show(inputArray,outputSequence,Math.min(inputArray.getRank()-1,2),inputSequence.getName()+" shifted", isHeadLess() );

        }else{
            inputArray =  ArrayUtils.roll(inputArray).copy();
            IcyImager.show(inputArray,outputSequence,inputSequence.getName()+" shifted", isHeadLess() );

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
