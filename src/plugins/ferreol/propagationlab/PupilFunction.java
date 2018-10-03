/**
 *
 */
package plugins.ferreol.propagationlab;

import javax.swing.JSeparator;

import icy.sequence.MetaDataUtil;
import icy.sequence.Sequence;
import icy.util.OMEUtil;
import loci.formats.ome.OMEXMLMetadataImpl;
import microTiPi.microUtils.Zernike;
import mitiv.array.ArrayFactory;
import mitiv.array.ArrayUtils;
import mitiv.array.Double3D;
import mitiv.base.Shape;
import mitiv.utils.MathUtils;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzLabel;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.ferreol.demics.ToolTipText;
import plugins.mitiv.io.IcyImager;

/**
 * @author ferreol
 *
 */
public class PupilFunction extends EzPlug implements Block, EzStoppable {
    protected EzVarInteger     npix;     // size of the output (npix x npix)

    // optical parameters
    protected EzVarDouble     dxy_nm;  //  pixels size in (x,y) in nm
    protected EzVarDouble     lambda;         //  wavelength
    protected EzVarDouble     ni;             //  refractive index of the immersion index
    protected EzVarDouble     na;             //  numerical aperture


    protected EzVarDouble defocus;
    protected EzVarDouble astigmatism0;
    protected EzVarDouble astigmatism45;
    protected EzVarDouble comaX;
    protected EzVarDouble comaY;
    protected EzVarDouble spherical;

    protected EzVarBoolean    fftshift;



    protected Shape pupilShape;

    private EzVarSequence   pupilout;
    int Nx, Ny;

    protected double[] zernike; // Zernike polynomials basis


    protected Sequence pupilSequence;
    @Override
    protected void initialize() {
        npix = new EzVarInteger("num pixel", 128,1, Integer.MAX_VALUE ,1);
        dxy_nm = new EzVarDouble("dxy(nm):",64.5,0., Double.MAX_VALUE,1.);
        lambda = new EzVarDouble( "\u03BB(nm):",540.,10.,15000.,5);
        ni = new EzVarDouble("ni:",1.,1.,2.,0.1);
        na= new EzVarDouble("na:",1.,0.,2.,0.1);
        fftshift = new EzVarBoolean("FFT shift", false);

        defocus = new EzVarDouble("defocus",0,-2, 2,.1);
        astigmatism0 = new EzVarDouble("astigmatism 0°",0,-2, 2,.1);
        astigmatism45 = new EzVarDouble("astigmatism 45°",0,-2, 2,.1);
        comaX = new EzVarDouble("vertical coma",0,-2, 2,.1);
        comaY = new EzVarDouble("horizontal coma",0,-2, 2,.1);
        spherical = new EzVarDouble("spherical",0,-2, 2,.1);


        addEzComponent(npix);
        npix.setToolTipText("number of pixels along each dimension");
        addEzComponent(dxy_nm);
        dxy_nm.setToolTipText(ToolTipText.doubleDxy);
        addEzComponent(lambda);
        lambda.setToolTipText(ToolTipText.doubleLambda);
        addEzComponent(ni);
        ni.setToolTipText(ToolTipText.doubleNi);
        addEzComponent(na);
        na.setToolTipText(ToolTipText.doubleNa);
        addEzComponent(fftshift);
        fftshift.setToolTipText("Swap quadrant to center the 0 frequency");


        addComponent(new JSeparator(JSeparator.VERTICAL));
        addEzComponent(new EzLabel("Aberration"));
        addEzComponent(defocus);
        addEzComponent(astigmatism0);
        addEzComponent(astigmatism45);
        addEzComponent(comaX);
        addEzComponent(comaY);
        addEzComponent(spherical);

        if (isHeadLess()) {
            pupilout = new EzVarSequence("generated plane wave");
        }
    }

    @Override
    protected void execute() {
        Nx = Ny =  npix.getValue(true);
        pupilShape = new Shape(2,Nx, Ny);
        Double3D   pupilArray =  Double3D.create(pupilShape);

        double[]  aberrationCoef =  new double[11];
        aberrationCoef[3] = defocus.getValue();
        aberrationCoef[4] = astigmatism45.getValue();
        aberrationCoef[5] = astigmatism0.getValue();
        aberrationCoef[6] = comaX.getValue();
        aberrationCoef[7] = comaY.getValue();
        aberrationCoef[10] = spherical.getValue();


        double pixsize = dxy_nm.getValue()*1E-9 ; //en metres
        double wavelenght = lambda.getValue() *1E-9;// en metres
        double k = 2.*Math.PI/wavelenght*ni.getValue();
        double radius =  (pixsize*Nx*na.getValue())/wavelenght;



        zernike = Zernike.zernikeArray(11, Nx, Ny, radius, true,false);
        zernike = MathUtils.gram_schmidt_orthonormalization(zernike, Nx, Ny, 11);
        Double3D zernikeArray = ArrayFactory.wrap(zernike, Nx, Ny, 11);


        for (int nx = 0; nx < Nx; nx++) {
            for (int ny = 0; ny < Ny; ny++) {
                double modulus = zernikeArray.get(nx, ny, 0);
                if (modulus!=0){
                    pupilArray.set(0, nx, ny, modulus);

                    double phase = 0;
                    for(int nz=3;nz<11;nz++){
                        phase += aberrationCoef[nz] *zernikeArray.get(nx, ny, nz);
                    }
                    pupilArray.set(1,nx, ny, phase);
                }
            }
        }

        if (fftshift.getValue()){
            int[] off = new int[3];
            off[1] = -(pupilArray.getDimension(1)/2);
            off[2] = -(pupilArray.getDimension(2)/2);
            pupilArray = (Double3D) ArrayUtils.roll(pupilArray,off).copy();
        }

        pupilSequence  = new Sequence();

        ome.xml.meta.OMEXMLMetadata newMetdat = MetaDataUtil.createMetadata("plane wave");
        newMetdat.setPixelsPhysicalSizeX(OMEUtil.getLength(dxy_nm.getValue()*1E-3), 0);
        newMetdat.setPixelsPhysicalSizeY(OMEUtil.getLength(dxy_nm.getValue()*1E-3), 0);
        newMetdat.setLaserWavelength(OMEUtil.getLength(lambda.getValue()*1E-3), 0, 0);
        newMetdat.setObjectiveSettingsRefractiveIndex(ni.getValue(), 0 );
        pupilSequence.setMetaData((OMEXMLMetadataImpl) newMetdat); //FIXME may not working now

        //   IcyImager.show(ArrayFactory.wrap(zernike, Nx, Ny, 11),pupilSequence,"zernike" ,isHeadLess() );

        IcyImager.show(pupilArray,pupilSequence,0,"Pupil" ,isHeadLess() );
        //    pupilSequence.setChannelName(0, "Real part");
        //    pupilSequence.setChannelName(1, "Imaginary part");

        if (isHeadLess()) {
            pupilout.setValue(pupilSequence);
        }

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

    /* (non-Javadoc)
     * @see plugins.adufour.ezplug.EzPlug#clean()
     */
    @Override
    public void clean() {
        // TODO Auto-generated method stub

    }


}
