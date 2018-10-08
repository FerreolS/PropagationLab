/**
 *
 */
package plugins.ferreol.PropagationLab;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import icy.plugin.interface_.PluginBundled;
import icy.sequence.Sequence;
import mitiv.array.ArrayUtils;
import mitiv.array.Double1D;
import mitiv.array.Double2D;
import mitiv.array.Double3D;
import mitiv.array.DoubleArray;
import mitiv.base.Shape;
import plugins.adufour.blocks.lang.Block;
import plugins.adufour.blocks.util.VarList;
import plugins.adufour.ezplug.EzPlug;
import plugins.adufour.ezplug.EzStoppable;
import plugins.adufour.ezplug.EzVarBoolean;
import plugins.adufour.ezplug.EzVarDouble;
import plugins.adufour.ezplug.EzVarInteger;
import plugins.adufour.ezplug.EzVarSequence;
import plugins.adufour.ezplug.EzVarText;
import plugins.ferreol.demics.ToolTipText;
import plugins.mitiv.io.Icy2TiPi;
import plugins.mitiv.io.IcyImager;

/**
 * @author ferreol
 *
 */
public class Pupil2PSF extends EzPlug implements Block, EzStoppable, PluginBundled {
    protected EzVarInteger     ezNz;
    protected EzVarDouble     dxy_nm;  //  pixels size in (x,y) in nm
    protected EzVarDouble     dz_nm;  //  pixels size in z in nm

    protected EzVarDouble     lambda;         //  wavelength
    protected EzVarDouble     ni;             //  refractive index of the immersion index
    protected EzVarDouble     na;             //  numerical aperture
    private EzVarSequence     pupil;
    protected EzVarBoolean    fftshiftout;
    protected EzVarBoolean    fftshiftin;

    protected EzVarText       inputOption;  // Combobox
    protected final static String[] inputOptions = new String[]{"Cartesian","Polar"};



    private EzVarSequence   PSF;

    protected Sequence pupilSequence;
    protected Sequence PSFSequence;

    protected Double3D pupilArray;
    protected static final double DEUXPI = 2*Math.PI;

    protected int Nxy;
    protected int Nz;
    protected double[] psi; // defocus function
    protected double dxy;
    private Shape psfShape;
    private Double3D psfArray;
    private Double dz;


    protected Shape psf2DShape;

    @Override
    protected void initialize() {
        if (!isHeadLess()) {
            getUI().setParametersIOVisible(false);
            //    getUI().setActionPanelVisible(false);
        }
        inputOption = new EzVarText(      "pupil  type:", inputOptions, false);

        pupil = new EzVarSequence("Pupil function");
        fftshiftin = new EzVarBoolean("FFT shift intput", false);

        ezNz = new EzVarInteger("number of depth", 64,1, Integer.MAX_VALUE ,1);
        dz_nm = new EzVarDouble("dz (nm):",64.5,0., Double.MAX_VALUE,1.);
        dxy_nm = new EzVarDouble("dxy(nm):",64.5,0., Double.MAX_VALUE,1.);
        lambda = new EzVarDouble( "\u03BB(nm):",540.,10.,15000.,5);
        ni = new EzVarDouble("refractive index:",1.,1.,2.,0.1);
        na= new EzVarDouble("NA:",1.,0.,2.,0.1);
        fftshiftout = new EzVarBoolean("FFT shift output", false);

        addEzComponent(pupil);
        addEzComponent(inputOption);
        addEzComponent(fftshiftin);
        addEzComponent(ezNz);
        ezNz.setToolTipText("Number of depth plane");
        addEzComponent(dz_nm);
        dxy_nm.setToolTipText(ToolTipText.doubleDz);
        addEzComponent(dxy_nm);
        dxy_nm.setToolTipText(ToolTipText.doubleDxy);
        addEzComponent(lambda);
        lambda.setToolTipText(ToolTipText.doubleLambda);
        addEzComponent(ni);
        ni.setToolTipText(ToolTipText.doubleNi);
        addEzComponent(na);
        na.setToolTipText(ToolTipText.doubleNa);
        addEzComponent(fftshiftout);

        if (isHeadLess()) {
            PSF = new EzVarSequence("PSF");
        }


    }

    @Override
    protected void execute() {

        pupilSequence = pupil.getValue();
        if (pupilSequence==null){
            return;
        }

        PSFSequence = new Sequence();

        DoubleArray tmpArray = Icy2TiPi.sequenceToArray(pupilSequence).toDouble();

        if ((tmpArray.getRank()!=3)||((tmpArray.getRank()==3)&&(tmpArray.getDimension(2)!=2))){
            throw new IllegalArgumentException("Pupil must be complex 2D");
        }
        if (tmpArray.getDimension(0)!=tmpArray.getDimension(1)){
            throw new IllegalArgumentException("Pupil must be square");
        }

        pupilArray = (Double3D) tmpArray.movedims(2, 0).toDouble().copy();
        Nxy = pupilArray.getDimension(1);

        if(inputOption.getValue()==inputOptions[0]){
            cartesian2polar(pupilArray);
        }

        if(fftshiftin.getValue()){

            int[] off = new int[3];
            off[1] = -(pupilArray.getDimension(1)/2);
            off[2] = -(pupilArray.getDimension(2)/2);
            pupilArray = (Double3D) ArrayUtils.roll(pupilArray,off).copy();

        }

        dz = dz_nm.getValue()*1E-9;
        Nz = ezNz.getValue();
        psfShape = new Shape(Nxy,Nxy,Nz);
        psfArray = Double3D.create(psfShape);



        psi = new double[Nxy*Nxy];

        dxy = dxy_nm.getValue()*1E-9;

        psf2DShape = new Shape(Nxy, Nxy);

        computeDefocus();
        computePsf();

        if( fftshiftout.getValue()){
            IcyImager.show(ArrayUtils.roll(psfArray),PSFSequence,"PSF",isHeadLess());
        }else{
            IcyImager.show(psfArray,PSFSequence,"PSF",isHeadLess());
        }
        //      PSFSequence.getFirstViewer().getLut().getLutChannel(0).setColorMap(new IceColorMap(),false);

        if (isHeadLess()) {
            PSF.setValue(PSFSequence);
        }
    }

    /**
     * @param pupilArray2
     */
    private void cartesian2polar(DoubleArray inArray) {
        for(int i=0;i<inArray.getNumber();i=i+2){
            double re = ((Double1D)inArray.toDouble().as1D()).get(i);
            double im = ((Double1D)inArray.toDouble().as1D()).get(i+1);
            ((Double1D)inArray.toDouble().as1D()).set(i,Math.sqrt(re*re+im*im));
            ((Double1D)inArray.toDouble().as1D()).set(i+1,Math.atan2(im,re));
        }
    }

    /**
     * Compute the defocus aberration ψ of the phase pupil
     * <p>
     */
    public void computeDefocus()
    {
        double lambda_ni2 = Math.pow(ni.getValue()/(lambda.getValue()*1E-9),2);
        double scale_x = 1/(Nxy*dxy);
        double scale_y = 1/(Nxy*dxy);
        double q, rx, ry;
        for (int ny = 0; ny < Nxy; ny++)
        {
            if(ny > Nxy/2)
            {
                ry = Math.pow(scale_y*(ny - Nxy) , 2);
            }
            else
            {
                ry = Math.pow(scale_y*ny, 2);
            }

            for (int nx = 0; nx < Nxy; nx++)
            {
                int nxy = nx + ny*Nxy;
                if ((pupilArray.get(0,nx, ny)!=0))
                {
                    if(nx > Nxy/2)
                    {
                        rx = Math.pow(scale_x*(nx - Nxy), 2);
                    }
                    else
                    {
                        rx = Math.pow(scale_x*nx, 2);
                    }

                    q = lambda_ni2 - rx - ry;

                    if (q < 0.0)
                    {
                        psi[nxy] = 0;
                        pupilArray.set(0,nx, ny,0);
                        pupilArray.set(1,nx, ny,0);
                    }
                    else
                    {
                        psi[nxy] = Math.sqrt(q);
                    }
                }
            }
        }
    }
    public void computePsf(){

        psfArray = Double3D.create( psfShape);

        final double PSFnorm = 1.0/(Nxy*Nxy*Nz);
        final int Npix = Nxy*Nxy;

        // double[] flattenPupil = pupilArray.flatten();
        if(true){
            int threads = Runtime.getRuntime().availableProcessors();
            ExecutorService service = Executors.newFixedThreadPool(threads);

            List<Future<GetPsfParaOut>> futures = new ArrayList<Future<GetPsfParaOut>>();
            for ( int iz = 0; iz < Nz; iz++)
            {
                final int iz1 = iz;
                Callable<GetPsfParaOut> callable = new Callable<GetPsfParaOut>() {
                    @Override
                    public GetPsfParaOut call() throws Exception {
                        GetPsfParaOut output = new GetPsfParaOut(Npix,iz1);
                        double defoc_scale;
                        double phasePupil;
                        double[] A = new double[2*Npix];

                        if (iz1 > Nz/2)
                        {
                            defoc_scale = DEUXPI*(iz1 - Nz)*dz;
                        }
                        else
                        {
                            defoc_scale = DEUXPI*iz1*dz;
                        }

                        for (int in = 0; in < Npix; in++)
                        {
                            phasePupil = pupilArray.as1D().get(2*in+1) + defoc_scale*psi[in];
                            double mod = pupilArray.as1D().get(2*in);

                            A[2*in] = mod*Math.cos(phasePupil);
                            A[2*in + 1] =mod*Math.sin(phasePupil);
                        }
                        /* Fourier transform of the pupil function A(z) */

                        DoubleFFT_2D   FFT2D = new DoubleFFT_2D(Nxy, Nxy);

                        FFT2D.complexForward(A);

                        for (int in = 0; in < Npix; in++)// FIXME use ShapedArray
                        {
                            ((double[])output.outPsf)[in] = (A[2*in]*A[2*in] + A[2*in+1]*A[2*in+1])*PSFnorm ;
                        }
                        return output;
                    }
                };
                futures.add(service.submit(callable));
            }

            service.shutdown();

            for (Future<GetPsfParaOut> future : futures) {
                GetPsfParaOut output;
                try {
                    output = future.get();
                    psfArray.slice(output.idxz).assign(Double2D.wrap((double[])output.outPsf, psf2DShape));
                } catch (InterruptedException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                } catch (ExecutionException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }else{
            DoubleFFT_2D FFT2D = new DoubleFFT_2D(Nxy, Nxy);

            for ( int iz = 0; iz < Nz; iz++)
            {
                double defoc_scale;
                double phasePupil;
                double[] A = new double[2*Npix];

                if (iz > Nz/2)
                {
                    defoc_scale = DEUXPI*(iz - Nz)*dz;
                }
                else
                {
                    defoc_scale = DEUXPI*iz*dz;
                }
                for (int in = 0; in < Npix; in++)
                {
                    phasePupil = pupilArray.as1D().get(2*in+1) + defoc_scale*psi[in];
                    double mod = pupilArray.as1D().get(2*in);

                    A[2*in] = mod*Math.cos(phasePupil);
                    A[2*in + 1] =mod*Math.sin(phasePupil);
                }
                /* Fourier transform of the pupil function A(z) */
                FFT2D.complexForward(A);

                for (int iy = 0; iy < Nxy; iy++){
                    for (int ix = 0; ix < Nxy; ix++){
                        int in = (ix+Nxy*iy);

                        psfArray.set(ix, iy, iz, (A[2*in]*A[2*in] + A[2*in+1]*A[2*in+1])*PSFnorm);
                    }
                }
            }
        }
    }





    private class GetPsfParaOut{
        Object  outPsf;
        int idxz;
        public GetPsfParaOut(int nPix, int iz){
            idxz = iz;
            outPsf = new double[2*nPix];
        }
    }


    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareInput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareInput(VarList inputMap) {
        initialize();
        inputMap.add("pupil function", pupil.getVariable());
        inputMap.add("input type", inputOption.getVariable());
        inputMap.add("FFTshift intput", fftshiftin.getVariable());
        inputMap.add("number of depth",ezNz.getVariable());
        inputMap.add("dz (nm)", dz_nm.getVariable());
        inputMap.add("dxy (nm)", dxy_nm.getVariable());
        inputMap.add( "\u03BB(nm)",lambda.getVariable());
        inputMap.add("Refractive index",ni.getVariable());
        inputMap.add("NA",na.getVariable());
        inputMap.add("FFTshift output",fftshiftout.getVariable());



    }

    /* (non-Javadoc)
     * @see plugins.adufour.blocks.lang.Block#declareOutput(plugins.adufour.blocks.util.VarList)
     */
    @Override
    public void declareOutput(VarList outputMap) {
        outputMap.add("PSF", PSF.getVariable());
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
