package plugins.ferreol.propagationlab;

import icy.sequence.Sequence;
import mitiv.array.Double3D;
import mitiv.base.Shape;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.ferreol.demics.ToolTipText;
import plugins.mitiv.io.IcyImager;

public class PlaneWave extends EzPlug implements Block {

    protected EzVarInteger     npix;     // size of the output (npix x npix)

    // optical parameters
    protected EzVarDouble     dxy_nm;  //  pixels size in (x,y) in nm
    protected EzVarDouble     lambda;         //  wavelength
    protected EzVarDouble     ni;             //  refractive index of the immersion index
    protected EzVarDouble     anglex , angley;

    protected Shape outputShape;

    private EzVarSequence   outputWave=null;

    protected Sequence outputSequence;

    int Nx, Ny;

    @Override
    protected void initialize() {
        // TODO Auto-generated by Icy4Eclipse

        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        npix = new EzVarInteger("num pixel", 128,1, Integer.MAX_VALUE ,1);
        dxy_nm = new EzVarDouble("dxy(nm):",64.5,0., Double.MAX_VALUE,1.);
        lambda = new EzVarDouble( "\u03BB(nm):",540.,10.,15000.,5);
        ni = new EzVarDouble("ni:",1.,1.,2.,0.1);
        anglex = new EzVarDouble("angle x:",0.,-90,90,5);
        angley = new EzVarDouble("angle y:",0.,-90,90,5);

        addEzComponent(npix);
        npix.setToolTipText("number of pixels along each dimension");
        addEzComponent(dxy_nm);
        dxy_nm.setToolTipText(ToolTipText.doubleDxy);
        addEzComponent(lambda);
        lambda.setToolTipText(ToolTipText.doubleLambda);
        addEzComponent(ni);
        ni.setToolTipText(ToolTipText.doubleNi);
        addEzComponent(anglex);
        anglex.setToolTipText("incidence angle along the first dimension");
        addEzComponent(angley);
        angley.setToolTipText("incidence angle along the first dimension");

    }

    @Override
    protected void execute() {
        // TODO Auto-generated by Icy4Eclipse
        Nx = Ny =  npix.getValue(true);
        outputShape = new Shape(2,Nx, Ny);
        Double3D   imgArray =  Double3D.create(outputShape);
        for (int nx = 0; nx < Nx; nx++) {
            for (int ny = 0; ny < Ny; ny++) {
                double phase = 2.*Math.PI / lambda.getValue()*dxy_nm.getValue()*(nx*Math.sin(anglex.getValue()/180*Math.PI)+ny*Math.sin(angley.getValue()/180*Math.PI));
                imgArray.set(0, nx, ny, Math.cos(phase));
                imgArray.set(1, nx, ny, Math.sin(phase));
            }
        }
        IcyImager.show(imgArray,outputSequence,0,"*" ,isHeadLess() );
        // MessageDialog.showDialog("PlaneWave is working fine !");
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
        initialize();
        inputMap.add("dxy", dxy_nm.getVariable());
        inputMap.add("npix", npix.getVariable());
        inputMap.add("ni", ni.getVariable());
        inputMap.add("lambda", lambda.getVariable());
        inputMap.add("angle X", anglex.getVariable());
        inputMap.add("angle Y", angley.getVariable());



    }

    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        // TODO Auto-generated method stub
        outputMap.add("outputWave", outputWave.getVariable());
    }
}
