package coreservlets;


import java.io.BufferedReader;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.Serializable;
import java.io.StringReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import javax.faces.application.FacesMessage;
import javax.faces.bean.ManagedBean;
import javax.faces.bean.SessionScoped;
import javax.faces.context.FacesContext;
import javax.servlet.ServletContext;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;
import javax.swing.JOptionPane;

import org.apache.poi.hssf.usermodel.HSSFRow;
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.primefaces.context.RequestContext;
import org.primefaces.model.DefaultStreamedContent;
import org.primefaces.model.StreamedContent;
import org.primefaces.model.chart.CartesianChartModel;
import org.primefaces.model.chart.LineChartSeries;





@ManagedBean
@SessionScoped

public class TuScatBean  implements Serializable{
	
	private String message;
	private CartesianChartModel linearModel=new CartesianChartModel(); 
	private CartesianChartModel linearModel2=new CartesianChartModel(); 
	private CartesianChartModel linearModel3=new CartesianChartModel(); 
	private CartesianChartModel linearModel4=new CartesianChartModel(); 
	private CartesianChartModel linearModel1=new CartesianChartModel(); 
	
	private HashMap<String, String> s11min=new HashMap<String, String>();
	
	
	private String rrefractiveIndexTextBox;
	private String irefractiveIndexTextBox;
	private String waveLengthText;
	
	private String selectcomboTypesShap;
	private String selectParameter;
	private String selectparameter2;
	
	
	
	private String combosizeDis;
	private String selectcombosizeDis;
	private String radiusTextBox;
	private String lowestPartRadiusText;
	private String highestPartRadiusText;
	private String sigmaText;
	private String modalradiusText;
	
	private String extinctionCoefficentValue;
	private String scattereingCoefficentValue;
	private String absorptionCoefficentValue;
	private String singleScatteringValue;
	private String assymmetricLabeltValue;
	private boolean boolsphere;
	private boolean boolshape;
	private boolean boolshapenot;
	private String labelalfasigma="Alfa";
	
	
	
	static final double PI =3.14159265e+0;
	static Complex refrelk,refmed;
	static Complex s1k[]=new Complex[200];
	static Complex s2k[]=new Complex[200];
	static double xk;
	static int nangk=91,nan,j;
	static double refre,refim,s11,s12,pol,s33,s34,ang,dang,s11nor,rad,wavel;
	
	
	/*Data Arrays*/
	static double angle[]=new double[181];
	static double ss11[]=new double[181];
	static double ss12[]=new double[181];
	static double ss33[]=new double[181];
	static double ss34[]=new double[181];
	static double ss11Com[]=new double[181];
	static double ss12Com[]=new double[181];
	static double ss33Com[]=new double[181];
	static double ss34Com[]=new double[181];
	static double ss11Min[]=new double[181];
	static double ss12Min[]=new double[181];
	static double ss33Min[]=new double[181];
	static double ss34Min[]=new double[181];
	static double d[]=new double[181];
	static double qscatxt;
	static double qbacktxt;
	static double qabstxt;
	static double albedotxt;
	static double qexttxt;
	static double gtxt;
	static double qprtxt;
	static double gscatxt;
	//static double valueGet=0.0;
	static int np;
	static boolean expBool=false;
	static double errorMin=-1.0;
	
	
	
	//experimental
	
	public boolean sizeRadioButton;
	
	public boolean refractiveRadioButton;
	
	public String fixedrrefractiveIndexTextBox;
	public String fixedirefractiveIndexTextBox;
	
	public boolean realindexRadioButton;
	public boolean imagindexRadioButton;
	
	public boolean sphereRadioButton;
	public boolean spheroidRadioButton;
	public boolean cyclinderRadioButton;
	
	
	public String selectsizeRef;
	
	public boolean boolsizeRI;
	public boolean boolsizeRISpCy;
	public boolean boolsizeRISproid;
	public boolean boolsizeABCD;
	public boolean boolsizeABCD1;
	public boolean boolsizeABCDspheroidab;
	public boolean boolsizeABCDspheroiddl;
	public boolean boolrefrctiveSphere;
	public boolean boolrefrctiveCylinder;
	public boolean boolrefrctiveSpheroid;
	public boolean boolrefractive;
	public boolean boolrefPara;
	
	public String minRadiusTextBox;
	public String maxRadiusTextBox;
	public String stepRadiusTextBox;
	
	
	private String accuracyTextBox;
	//private String a_bOrC_LTextBox;
	private String selectcomboTypesPara;
	
	private String eqvRadiusTextBox;
	private String minAccuracyTextBox;
	private String maxAccuracyTextBox;
	private String stepAccuracyTextBox;
	
	private String displayTextABVE;
	private String displayTextComBo;
	private String displayTextComBoEqui;
	private String displayTextRealImag;
	private String displayTextLegen;
	private String displayTextLast;
	
	private String displayexpectedRange;
	
	private String selectcomboTypesShapImg;
	
	
	private String  radiusImgRealTextBox;
	private String  minirefractiveIndexTextBox;
	private String  maxrefractiveIndexTextBox;
	private String  steprefractiveIndexTextBox;
	private String a_bOrC_LTextBox;
	
	private StreamedContent file;
	public FacesMessage massage;
	public String userid;
	
	private List<TuScatPojo> listofdata=new ArrayList<TuScatPojo>();
	
	static String minimunError;
	
	/*public TuScatBean()
	{
			FacesContext context = FacesContext.getCurrentInstance();
			HttpServletRequest reques = (HttpServletRequest)context.getExternalContext().getRequest();  
			HttpSession session = reques.getSession(true);
			userid = (String) session.getAttribute("emailid");
			
			FacesMessage massage = null;
			massage = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "Successfull login"+"..."+userid);
			FacesContext.getCurrentInstance().addMessage(null, massage);
			
				
	}*/
	
	
	
	////////////////////////////////////For Calculation of minimum Error//////////////////////////////////////
	public void minimumCalculationErrorMessage() throws NullPointerException
	{
		try
		{
			StringReader sr= new StringReader(minimunError);
			BufferedReader br = new BufferedReader(sr);
			ReadWithScanner readWithScanner=new ReadWithScanner();;
			String strLine;
			//Read File Line By Line
			String test="";
			while ((strLine = br.readLine()) != null)   {
			  // Print the content on the console
			  //System.out.println (strLine);
				test+=strLine+"\n";
				readWithScanner.processLine(strLine);
			}
			double radlk=0.0;
			double radhk=0.0;
			double stepk;
			int ntotk;
			//boolCalculation=true;
				if(sizeRadioButton==true)
					{
						if(boolsizeRI==true)
						{
							errorMin=-1.0;
							String strrefre=fixedrrefractiveIndexTextBox;
							refre=Double.parseDouble(strrefre);
							String strrefim=fixedirefractiveIndexTextBox;
							refim=Double.parseDouble(strrefim);
							String strrefmedreal="1";
							double real=Double.parseDouble(strrefmedreal);
							refmed=new Complex(real,0.0);
							refrelk=new Complex(refre,refim);
							double valueGetsizesp=0.0;
							ntotk=1; 
							stepk=0.05;
							String strxk1=waveLengthText;
							double lam=Double.parseDouble(strxk1);
							double minRadius= Double.parseDouble(minRadiusTextBox);
							double maxRadius= Double.parseDouble(maxRadiusTextBox);
							double stepRadius= Double.parseDouble(stepRadiusTextBox);
							for(double j=minRadius;j<=maxRadius;j=j+stepRadius)
							{
								xk=(2.0*PI*j*real)/Double.parseDouble(strxk1);
								main_Function();
								double sum=0;
								for(int i=10;i<171;i++)
								{
									String seq=""+i;
									String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
									d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
									sum=sum+d[i];
									System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
								}
								
								double rmserror;
								rmserror=Math.sqrt(sum/161);
								if(errorMin==-1.0)
								{
									errorMin=rmserror;
									valueGetsizesp=j;
								}
								else if(rmserror<errorMin)
								{
									errorMin=rmserror;
									valueGetsizesp=j;
								}
							}
							
							RequestContext context = RequestContext.getCurrentInstance();  
							FacesMessage msg1 = null;
							msg1 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetsizesp);
							FacesContext.getCurrentInstance().addMessage(null, msg1);
							//System.out.println("Refractive Index(Real) :"+fixedrrefractiveIndexTextBox+"\n");
							//System.out.println("Refractive Index(Imaginary) :"+fixedirefractiveIndexTextBox+"\n");
							//System.out.println("Wave Length Value :"+waveLengthText+"\n");
							//System.out.println("Radius :"+valueGet+"\n");
			
						}/////////////////// End of if for D/L Ratio
								
						
						
								///////////////Start of Cylinder
						else if(boolsizeRISpCy==true)
						{
							
							////D/L
								if(boolsizeABCD==true)
									{
									System.out.println("Hiiiiiiiiiiiiiiiii this is cylinder with D/L ration");
									//System.out.println("Hiiiiiiiiiiiiiiiii");
										errorMin=-1.0;
										String strrefre=fixedrrefractiveIndexTextBox;
										refre=Double.parseDouble(strrefre);
										String strrefim=fixedirefractiveIndexTextBox;
										refim=Double.parseDouble(strrefim);
										String strrefmedreal="1";
										double real=Double.parseDouble(strrefmedreal);
										refmed=new Complex(real,0.0);
										refrelk=new Complex(refre,refim);
										ntotk=1; 
										stepk=0.05;
										String strxk1=waveLengthText;
										double lam=Double.parseDouble(strxk1);
										double axi=Double.parseDouble(eqvRadiusTextBox);
										double minRadius= Double.parseDouble(minAccuracyTextBox);
										double maxRadius= Double.parseDouble(maxAccuracyTextBox);
										double stepRadius= Double.parseDouble(stepAccuracyTextBox);
										String strddelt=accuracyTextBox;
										double ddelt=Double.parseDouble(strddelt);
										double valueGetsizecydl=0.0;
										System.out.println("Hi init");
										for(double j=minRadius;j<=maxRadius;j=j+stepRadius)
										{
								
											radlk=0.9999999*axi;
											radhk=1.0000001*axi;
											String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,j,ddelt,-1,-2,lam);//eps a/b or d/l ratio
											ss11=nonSphericalClass.ss11;
											ss12=nonSphericalClass.ss12;
											ss33=nonSphericalClass.ss33;
											ss34=nonSphericalClass.ss34;
											qexttxt=nonSphericalClass.cextin;
											qscatxt=nonSphericalClass.cscat;
											qabstxt=nonSphericalClass.cabsin;
											albedotxt=nonSphericalClass.walb;
											gtxt=nonSphericalClass.asymm;
											double sum=0;
											System.out.println("Hi Claculating");
											for(int i=10;i<171;i++)
											{
												
												String seq=""+i;
												String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
												d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
												sum=sum+d[i];
												System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
											double rmserror;
											rmserror=Math.sqrt(sum/181);
											if(errorMin==-1.0)
											{
												errorMin=rmserror;
												valueGetsizecydl=j;
											}
											else if(rmserror<errorMin)
											{
												errorMin=rmserror;
												valueGetsizecydl=j;
							
											}
										}
										RequestContext context2 = RequestContext.getCurrentInstance();  
										FacesMessage msg2 = null;
										msg2 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetsizecydl);
										FacesContext.getCurrentInstance().addMessage(null, msg2);
							
										//System.out.println(""+"Wave Length Value :"+waveLengthText+"\n");
										//System.out.println(""+"Refractive Index(Real) :"+fixedrrefractiveIndexTextBox+"\n");
										//System.out.println(""+"Refractive Index(Imaginary) :"+fixedirefractiveIndexTextBox+"\n");
										//System.out.println(""+"Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
										//System.out.println(""+"Equivalent Sphere Radius :"+eqvRadiusTextBox+"\n");
										//System.out.println(""+"A/B:"+valueget+"\n");
							
									}/////////////////// End of if for D/L Ratio
							
							
							///Equivalent
								else if(boolsizeABCD1==true)
									{
									System.out.println("Hi This is Equivalent radious for Cylinder"+boolsizeRISpCy+""+boolsizeABCD1);
										errorMin=-1.0;
										String strrefre=fixedrrefractiveIndexTextBox;
										refre=Double.parseDouble(strrefre);
										String strrefim=fixedirefractiveIndexTextBox;
										refim=Double.parseDouble(strrefim);
										String strrefmedreal="1";
										double real=Double.parseDouble(strrefmedreal);
										refmed=new Complex(real,0.0);
										refrelk=new Complex(refre,refim);
										ntotk=1; 
										stepk=0.05;
										String strxk1=waveLengthText;
										double lam=Double.parseDouble(strxk1);
										double eps=Double.parseDouble(eqvRadiusTextBox);
										double minRadius= Double.parseDouble(minAccuracyTextBox);
										double maxRadius= Double.parseDouble(maxAccuracyTextBox);
										double stepRadius= Double.parseDouble(stepAccuracyTextBox);;
										String strddelt=accuracyTextBox;
										double ddelt=Double.parseDouble(strddelt);
										double valueGetsizecyeq=0.0;
										
										for(double j=minRadius;j<=maxRadius;j=j+stepRadius)
											{
											
											radlk=0.9999999*j;
											radhk=1.0000001*j;
											String ret=nonSphericalClass.nonsp(refre,refim,1,j,radlk,radhk,7,eps,ddelt,-1,-2,lam);
											ss11=nonSphericalClass.ss11;
											ss12=nonSphericalClass.ss12;
											ss33=nonSphericalClass.ss33;
											ss34=nonSphericalClass.ss34;
											qexttxt=nonSphericalClass.cextin;
											qscatxt=nonSphericalClass.cscat;
											qabstxt=nonSphericalClass.cabsin;
											albedotxt=nonSphericalClass.walb;
											gtxt=nonSphericalClass.asymm;
											double sum=0;
											for(int i=10;i<171;i++)
											{
												String seq=""+i;
												String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
												d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
												sum=sum+d[i];
												System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
											System.out.println("CylinEqr");
											double rmserror;
											rmserror=Math.sqrt(sum/181);
											if(errorMin==-1.0)
											{
												errorMin=rmserror;
												valueGetsizecyeq=errorMin;
											}
											else if(rmserror<errorMin)
											{
												errorMin=rmserror;
												valueGetsizecyeq=errorMin;
						
											}
											}
										RequestContext context = RequestContext.getCurrentInstance();  
										FacesMessage msg3 = null;
										msg3 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetsizecyeq);
										FacesContext.getCurrentInstance().addMessage(null, msg3);
							
										//System.out.println(""+"Wave Length Value :"+waveLengthText+"\n");
										//System.out.println(""+"Refractive Index(Real) :"+fixedrrefractiveIndexTextBox+"\n");
										//System.out.println(""+"Refractive Index(Imaginary) :"+fixedirefractiveIndexTextBox+"\n");
										//System.out.println(""+"Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
										//System.out.println(""+"A/B :"+eqvRadiusTextBox+"\n");
										//System.out.println(""+"Equivalent Sphere Radius:"+valueGet+"\n");
							
									}
								//// End of Else if for Equivalent Radius
						}
						//// End of Else if for boolsizeRISpCy oR Cylinder
						
						//////////////////////// Start of Spheroid
						else if(boolsizeRISproid==true)
						{
							
							System.out.println("boolsizeRISproid"+boolsizeRISproid);
								/////Start of A/B Ratio spheroid
								if(boolsizeABCDspheroidab==true)/////Start of A/B Ratio spheroid
								{
									System.out.println("boolsizeABCDspheroidab"+boolsizeABCDspheroidab);
									errorMin=-1.0;
									String strrefre=fixedrrefractiveIndexTextBox;
									refre=Double.parseDouble(strrefre);
									String strrefim=fixedirefractiveIndexTextBox;
									refim=Double.parseDouble(strrefim);
									String strrefmedreal="1";
									double real=Double.parseDouble(strrefmedreal);
									refmed=new Complex(real,0.0);
									refrelk=new Complex(refre,refim);
									ntotk=1; 
									stepk=0.05;
									String strxk1=waveLengthText;
									double lam=Double.parseDouble(strxk1);
									double axi=Double.parseDouble(eqvRadiusTextBox);
									double minRadius= Double.parseDouble(minAccuracyTextBox);
									double maxRadius= Double.parseDouble(maxAccuracyTextBox);
									double stepRadius= Double.parseDouble(stepAccuracyTextBox);
									String strddelt=accuracyTextBox;
									double ddelt=Double.parseDouble(strddelt);
									double valueGetsizeSpab=0.0;
									
									System.out.println(""+"Wave Length Value :"+lam+"\n");
									System.out.println(""+"axi :"+axi+"\n");
									System.out.println(""+"strrefre :"+strrefre+"\n");
									System.out.println(""+"strrefim :"+strrefim+"\n");
									System.out.println(""+"minRadius :"+minRadius+"\n");
									System.out.println(""+"maxRadius :"+maxRadius+"\n");
									System.out.println(""+"stepRadius :"+stepRadius+"\n");
									System.out.println(""+"ddelt :"+ddelt+"\n");
									
									
									
									
									
									for(double j=minRadius;j<=maxRadius;j=j+stepRadius)
									{
										radlk=0.9999999*axi;
										radhk=1.0000001*axi;
										String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,j,ddelt,-1,-1,lam);
								
										ss11=nonSphericalClass.ss11;
										ss12=nonSphericalClass.ss12;
										ss33=nonSphericalClass.ss33;
										ss34=nonSphericalClass.ss34;
										qexttxt=nonSphericalClass.cextin;
										qscatxt=nonSphericalClass.cscat;
										qabstxt=nonSphericalClass.cabsin;
										albedotxt=nonSphericalClass.walb;
										gtxt=nonSphericalClass.asymm;
										double sum=0;
										for(int i=10;i<171;i++)
											{
											String seq=""+i;
											String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
											d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
											sum=sum+d[i];
											System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
										double rmserror;
										rmserror=Math.sqrt(sum/181);
										if(errorMin==-1.0)
											{
											errorMin=rmserror;
											valueGetsizeSpab=j;
											}
										else if(rmserror<errorMin)
										{
											errorMin=rmserror;
											valueGetsizeSpab=j;
											
										}
									}
									RequestContext context = RequestContext.getCurrentInstance();  
									FacesMessage msg4 = null;
									msg4 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetsizeSpab);
									FacesContext.getCurrentInstance().addMessage(null, msg4);
						
									//System.out.println(""+"Wave Length Value :"+waveLengthText+"\n");
									//System.out.println(""+"Refractive Index(Real) :"+fixedrrefractiveIndexTextBox+"\n");
									//System.out.println(""+"Refractive Index(Imaginary) :"+fixedirefractiveIndexTextBox+"\n");
									//System.out.println(""+"Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
									//System.out.println(""+"Equivalent Sphere Radius :"+eqvRadiusTextBox+"\n");
									//System.out.println(""+"A/B:"+valueget+"\n");
						
								}
								//////////////End of A/B Ratio Spheroid
								/////////////Start of Eqv Radius spheroid
								else if(boolsizeABCDspheroiddl==true)
									{
									System.out.println("boolsizeRISproid Equivalent Radius"+boolsizeABCDspheroiddl);
									errorMin=-1.0;
									String strrefre=fixedrrefractiveIndexTextBox;
									refre=Double.parseDouble(strrefre);
									String strrefim=fixedirefractiveIndexTextBox;
									refim=Double.parseDouble(strrefim);
									String strrefmedreal="1";
									double real=Double.parseDouble(strrefmedreal);
									refmed=new Complex(real,0.0);
									refrelk=new Complex(refre,refim);
									ntotk=1; 
									stepk=0.05;
									String strxk1=waveLengthText;
									double lam=Double.parseDouble(strxk1);
									double axi=Double.parseDouble(eqvRadiusTextBox);
									double minRadius= Double.parseDouble(minAccuracyTextBox);
									double maxRadius= Double.parseDouble(maxAccuracyTextBox);
									double stepRadius= Double.parseDouble(stepAccuracyTextBox);
									String strddelt=accuracyTextBox;
									double ddelt=Double.parseDouble(strddelt);
									double valueGetsizeSpeq=0.0;
										for(double j=minRadius;j<=maxRadius;j=j+stepRadius)
										{
											radlk=0.9999999*axi;
											radhk=1.0000001*axi;
											String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,j,ddelt,-1,-1,lam);
											ss11=nonSphericalClass.ss11;
											ss12=nonSphericalClass.ss12;
											ss33=nonSphericalClass.ss33;
											ss34=nonSphericalClass.ss34;
											qexttxt=nonSphericalClass.cextin;
											qscatxt=nonSphericalClass.cscat;
											qabstxt=nonSphericalClass.cabsin;
											albedotxt=nonSphericalClass.walb;
											gtxt=nonSphericalClass.asymm;
											double sum=0;
											for(int i=10;i<171;i++)
											{
												String seq=""+i;
												String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
												d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
												sum=sum+d[i];
												System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
											double rmserror;
											rmserror=Math.sqrt(sum/181);
											if(errorMin==-1.0)
											{
												errorMin=rmserror;
												valueGetsizeSpeq=j;
											}
											else if(rmserror<errorMin)
											{
												errorMin=rmserror;
												valueGetsizeSpeq=j;
											}
										}
										RequestContext context = RequestContext.getCurrentInstance();  
										FacesMessage msg5 = null;
										msg5 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetsizeSpeq);
										FacesContext.getCurrentInstance().addMessage(null, msg5);
										
										//System.out.println(""+"Wave Length Value :"+waveLengthText+"\n");
										//System.out.println(""+"Refractive Index(Real) :"+fixedrrefractiveIndexTextBox+"\n");
										//System.out.println(""+"Refractive Index(Imaginary) :"+fixedirefractiveIndexTextBox+"\n");
										//System.out.println(""+"Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
										//System.out.println(""+"D/L :"+eqvRadiusTextBox+"\n");
										//System.out.println(""+"Equivalent Sphere Radius:"+valueGet+"\n");
									}    /////////////End of Eqv Radius spheroid
						}  ////////////////End Of Spheroid
						
					}//////////////End of Size

				
				////////////////////Start Of Refractive Index
				else if(refractiveRadioButton==true)
				{
					System.out.println("refractiveRadioButton"+refractiveRadioButton);/////////////Refractive Index and real part
					if(realindexRadioButton==true)
								{
						System.out.println("realindexRadioButton"+realindexRadioButton);/////////////Refractive Index and real part and Sphere
									if(boolrefrctiveSphere==true)
									{
										System.out.println("boolrefrctiveSphere"+boolrefrctiveSphere);
										errorMin=-1.0;
										double valueGetrefrealshp=0.0;
										double strrefremax=Double.parseDouble(maxrefractiveIndexTextBox);
										double strrefremin=Double.parseDouble(minirefractiveIndexTextBox);
										double strrefrestep=Double.parseDouble(steprefractiveIndexTextBox);
										for(double j=strrefremin;j<=strrefremax;j=j+strrefrestep)
										{
											refre=j;
											String strrefim=radiusImgRealTextBox;
											refim=Double.parseDouble(strrefim);
											double real=1.0;
											
											refmed=new Complex(real,0.0);
											refrelk=new Complex(refre,refim);
											ntotk=1; 
											stepk=0.05;
											String strxk1=waveLengthText;
											double lam=Double.parseDouble(strxk1);
											String strxk2 =radiusTextBox;
											xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
											main_Function();
											double sum=0;
											for(int i=10;i<171;i++)
											{
												String seq=""+i;
												String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
												d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
												sum=sum+d[i];
												System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
											double rmserror;
											rmserror=Math.sqrt(sum/181);
											if(errorMin==-1.0)
											{
												errorMin=rmserror;
												valueGetrefrealshp=j;
											}
											else if(rmserror<errorMin)
											{
												errorMin=rmserror;
												valueGetrefrealshp=j;
											}
										}
				
										RequestContext context = RequestContext.getCurrentInstance();  
										FacesMessage msg6 = null;
										msg6 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetrefrealshp);
										FacesContext.getCurrentInstance().addMessage(null, msg6);
										System.out.println("Wave Length Value :"+waveLengthText+"\n");
										System.out.println("Refractive Index(Real) :"+valueGetrefrealshp+"\n");
										System.out.println("Refractive Index(Imaginary) :"+radiusImgRealTextBox+"\n");
										System.out.println("Radius::"+radiusTextBox+"\n");
			
									}
									////////End of Sphere
						/////////////Refractive Index and real part and Chylinder
									else if(boolrefrctiveCylinder==true)
									{
										System.out.println("boolrefrctiveCylinder"+boolrefrctiveCylinder);
										double valueGetrefrealcyl=0.0;
										errorMin=-1.0;
										double strrefremax=Double.parseDouble(maxrefractiveIndexTextBox);
										double strrefremin=Double.parseDouble(minirefractiveIndexTextBox);
										double strrefrestep=Double.parseDouble(steprefractiveIndexTextBox);
										for(double j=strrefremin;j<=strrefremax;j=j+strrefrestep)
										{
											refre=j;
											String strrefim=radiusImgRealTextBox;
											refim=Double.parseDouble(strrefim);
											double real=1.0;
											refmed=new Complex(real,0.0);
											refrelk=new Complex(refre,refim);
											ntotk=1; 
											stepk=0.05;
											String strxk1=waveLengthText;
											double lam=Double.parseDouble(strxk1);
											String strxk2 =eqvRadiusTextBox;
											xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
											double axi=Double.parseDouble(strxk2);
											String streps=a_bOrC_LTextBox;
											double eps=Double.parseDouble(streps);
											String strddelt=accuracyTextBox;
											double ddelt=Double.parseDouble(strddelt);
											radlk=0.9999999*axi;
											radhk=1.0000001*axi;
											String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,eps,ddelt,-1,-2,lam);
										
											ss11=nonSphericalClass.ss11;
											ss12=nonSphericalClass.ss12;
											ss33=nonSphericalClass.ss33;
											ss34=nonSphericalClass.ss34;
											qexttxt=nonSphericalClass.cextin;
											qscatxt=nonSphericalClass.csca;
											qabstxt=nonSphericalClass.cabsin;
											albedotxt=nonSphericalClass.walb;
											gtxt=nonSphericalClass.asymm;
											double sum=0;
											for(int i=10;i<171;i++)
											{
												String seq=""+i;
												String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
												d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
												sum=sum+d[i];
												System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
											double rmserror;
											rmserror=Math.sqrt(sum/181);
											System.out.println("rmserror :"+rmserror+"\n");
											System.out.println("Real :"+j+"\n");
											
											
											if(errorMin==-1.0)
											{
												errorMin=rmserror;
												valueGetrefrealcyl=j;
											}
											else if(rmserror<errorMin)
											{
												errorMin=rmserror;
												valueGetrefrealcyl=j;
												
											}
										}
										RequestContext context = RequestContext.getCurrentInstance();  
										FacesMessage msg7 = null;
										msg7 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Real Estimate"+valueGetrefrealcyl);
										FacesContext.getCurrentInstance().addMessage(null, msg7);
				
										//System.out.println("Wave Length Value :"+waveLengthText+"\n");
										//System.out.println("Refractive Index(Real) :"+valueGetrefrealcyl+"\n");
										//System.out.println("Refractive Index(Imaginary) :"+radiusImgRealTextBox+"\n");
										//System.out.println("Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
										//System.out.println("D/L :"+a_bOrC_LTextBox+"\n");
										//System.out.println("Equivalent Sphere Radius:"+eqvRadiusTextBox+"\n");
				
									}
							////////End of Chylinder
						/////////////Start of Refractive Index and real part and Spheroid
									else if(boolrefrctiveSpheroid==true)
									{
										System.out.println("boolrefrctiveSpheroid"+boolrefrctiveSpheroid);
										double strrefremax=Double.parseDouble(maxrefractiveIndexTextBox);
										double strrefremin=Double.parseDouble(minirefractiveIndexTextBox);
										double strrefrestep=Double.parseDouble(steprefractiveIndexTextBox);
										double valueGetrefresphro=0.0;
										for(double j=strrefremin;j<=strrefremax;j=j+strrefrestep)
										{
											refre=j;
											errorMin=-1.0;
											String strrefim=radiusImgRealTextBox;
											refim=Double.parseDouble(strrefim);
											double real=1.0;
											refmed=new Complex(real,0.0);
											refrelk=new Complex(refre,refim);
											ntotk=1; 
											stepk=0.05;
											String strxk1=waveLengthText;
											double lam=Double.parseDouble(strxk1);
											String strxk2 =eqvRadiusTextBox;
											xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
											double axi=Double.parseDouble(strxk2);
											String streps=a_bOrC_LTextBox;
											double eps=Double.parseDouble(streps);
											String strddelt=accuracyTextBox;
											double ddelt=Double.parseDouble(strddelt);
											radlk=0.9999999*axi;
											radhk=1.0000001*axi;
											String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,eps,ddelt,-1,-1,lam);
						
											ss11=nonSphericalClass.ss11;
											ss12=nonSphericalClass.ss12;
											ss33=nonSphericalClass.ss33;
											ss34=nonSphericalClass.ss34;
											qexttxt=nonSphericalClass.cextin;
											qscatxt=nonSphericalClass.csca;
											qabstxt=nonSphericalClass.cabsin;
											albedotxt=nonSphericalClass.walb;
											gtxt=nonSphericalClass.asymm;
											double sum=0;
											for(int i=10;i<171;i++)
											{
												String seq=""+i;
												String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
												d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
												sum=sum+d[i];
												System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
											}
											double rmserror;
											rmserror=Math.sqrt(sum/181);
											if(errorMin==-1.0)
											{
												errorMin=rmserror;
												valueGetrefresphro=j;
											}
											else if(rmserror<errorMin)
											{
												errorMin=rmserror;
												valueGetrefresphro=j;
											}
										}
										RequestContext context = RequestContext.getCurrentInstance();  
										FacesMessage msg8 = null;
										msg8 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetrefresphro);
										FacesContext.getCurrentInstance().addMessage(null, msg8);
										System.out.println("Wave Length Value :"+waveLengthText+"\n");
										System.out.println("Refractive Index(Real) :"+valueGetrefresphro+"\n");
										System.out.println("Refractive Index(Imaginary) :"+radiusImgRealTextBox+"\n");
										System.out.println("Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
										System.out.println("A/B :"+a_bOrC_LTextBox+"\n");
										System.out.println("Equivalent Sphere Radius:"+eqvRadiusTextBox+"\n");
				
									}
						/////////////Start of Refractive Index and real part and Spheroid
								}	
					/////////////End of Refractive Index and real part
					
					
					
					
					
					
					
		/////////////Start of Refractive Index and Imaginary part
					else if(imagindexRadioButton==true)
					{
						System.out.println("imagindexRadioButton"+imagindexRadioButton);/////////////Start of Refractive Index and Imaginary part and sphere
							if(boolrefrctiveSphere==true)
							{
								System.out.println("boolrefrctiveSphere"+boolrefrctiveSphere);
								double valueGetrefimasph=0.0;
								errorMin=-1.0;
								double strrefremax=Double.parseDouble(maxrefractiveIndexTextBox);
								double strrefremin=Double.parseDouble(minirefractiveIndexTextBox);
								double strrefrestep=Double.parseDouble(steprefractiveIndexTextBox);
								for(double j=strrefremin;j<=strrefremax;j=j+strrefrestep)
								{
									String strrefre=radiusImgRealTextBox;
									refre=Double.parseDouble(strrefre);
									refim=j;
									double real=1.0;
									refmed=new Complex(real,0.0);
									refrelk=new Complex(refre,refim);
									ntotk=1; 
									stepk=0.05;
									String strxk1=waveLengthText;
									double lam=Double.parseDouble(strxk1);
									String strxk2 =radiusTextBox;
									xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
									main_Function();
									double sum=0;
									for(int i=10;i<171;i++)
									{
										String seq=""+i;
										String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
										d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
										sum=sum+d[i];
										System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
									}
									double rmserror;
									rmserror=Math.sqrt(sum/181);
									if(errorMin==-1.0)
									{
									errorMin=rmserror;
									valueGetrefimasph=j;
									}
									else if(rmserror<errorMin)
									{
									errorMin=rmserror;
									valueGetrefimasph=j;
									}
								}
		
							RequestContext context = RequestContext.getCurrentInstance();  
							FacesMessage msg9 = null;
							msg9 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetrefimasph);
							FacesContext.getCurrentInstance().addMessage(null, msg9);
							System.out.println("Wave Length Value :"+waveLengthText+"\n");
							System.out.println("Refractive Index(Real) :"+valueGetrefimasph+"\n");
							System.out.println("Refractive Index(Imaginary) :"+radiusImgRealTextBox+"\n");
							System.out.println("Radius::"+radiusTextBox+"\n");
		
						}
			/////////////End of Refractive Index and Imaginary part and sphere
			/////////////Start of Refractive Index and Imaginary part and cyclinder
						else if(boolrefrctiveCylinder==true)
						{
							System.out.println("boolrefrctiveCylinder"+boolrefrctiveCylinder);
							double strrefremax=Double.parseDouble(maxrefractiveIndexTextBox);
							double strrefremin=Double.parseDouble(minirefractiveIndexTextBox);
							double strrefrestep=Double.parseDouble(steprefractiveIndexTextBox);
							double valueGetrefimacyl=0.0;
							errorMin=-1.0;
							for(double j=strrefremin;j<=strrefremax;j=j+strrefrestep)
							{
								String strrefre=radiusImgRealTextBox;
								refre=Double.parseDouble(strrefre);
								refim=j;
								double real=1.0;
								refmed=new Complex(real,0.0);
								refrelk=new Complex(refre,refim);
								ntotk=1; 
								stepk=0.05;
								String strxk1=waveLengthText;
								double lam=Double.parseDouble(strxk1);
								String strxk2 =eqvRadiusTextBox;
								xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
								double axi=Double.parseDouble(strxk2);
								String streps=a_bOrC_LTextBox;
								double eps=Double.parseDouble(streps);
								String strddelt=accuracyTextBox;
								double ddelt=Double.parseDouble(strddelt);
								radlk=0.9999999*axi;
								radhk=1.0000001*axi;
								String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,eps,ddelt,-1,-1,lam);
		
								ss11=nonSphericalClass.ss11;
								ss12=nonSphericalClass.ss12;
								ss33=nonSphericalClass.ss33;
								ss34=nonSphericalClass.ss34;
								qexttxt=nonSphericalClass.cextin;
								qscatxt=nonSphericalClass.csca;
								qabstxt=nonSphericalClass.cabsin;
								albedotxt=nonSphericalClass.walb;
								gtxt=nonSphericalClass.asymm;
								double sum=0;
								for(int i=10;i<171;i++)
								{
									String seq=""+i;
									String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
									d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
									sum=sum+d[i];
									System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
								}
								double rmserror;
								rmserror=Math.sqrt(sum/181);
								if(errorMin==-1.0)
								{
									errorMin=rmserror;
									valueGetrefimacyl=j;
								}
								else if(rmserror<errorMin)
								{
									errorMin=rmserror;
									valueGetrefimacyl=j;
								}
							}
							RequestContext context = RequestContext.getCurrentInstance();  
							FacesMessage msg10 = null;
							msg10 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetrefimacyl);
							FacesContext.getCurrentInstance().addMessage(null, msg10);
		
							System.out.println("Wave Length Value :"+waveLengthText+"\n");
							System.out.println("Refractive Index(Real) :"+ radiusImgRealTextBox+"\n");
							System.out.println("Refractive Index(Imaginary) :"+ valueGetrefimacyl+"\n");
							System.out.println("Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
							System.out.println("A/B :"+a_bOrC_LTextBox+"\n");
							System.out.println("Equivalent Sphere Radius:"+eqvRadiusTextBox+"\n");
		
						}
			/////////////End of Refractive Index and Imaginary part and cyclinder
			/////////////Start of Refractive Index and Imaginary part and spheroid			
						else if(boolrefrctiveSpheroid==true)
						{
							System.out.println("boolrefrctiveSpheroid"+boolrefrctiveSpheroid);
							double strrefremax=Double.parseDouble(maxrefractiveIndexTextBox);
							double strrefremin=Double.parseDouble(minirefractiveIndexTextBox);
							double strrefrestep=Double.parseDouble(steprefractiveIndexTextBox);
							double valueGetrefimasphero=0.0;
							errorMin=-1.0;
							for(double j=strrefremin;j<=strrefremax;j=j+strrefrestep)
							{
								String strrefre=radiusImgRealTextBox;
								refre=Double.parseDouble(strrefre);
								refim=j;
								double real=1.0;
								refmed=new Complex(real,0.0);
								refrelk=new Complex(refre,refim);
								ntotk=1; 
								stepk=0.05;
								String strxk1=waveLengthText;
								double lam=Double.parseDouble(strxk1);
								String strxk2 =eqvRadiusTextBox;
								xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
								double axi=Double.parseDouble(strxk2);
								String streps=a_bOrC_LTextBox;
								double eps=Double.parseDouble(streps);
								String strddelt=accuracyTextBox;
								double ddelt=Double.parseDouble(strddelt);
								radlk=0.9999999*axi;
								radhk=1.0000001*axi;
								String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,eps,ddelt,-1,-1,lam);
		
								ss11=nonSphericalClass.ss11;
								ss12=nonSphericalClass.ss12;
								ss33=nonSphericalClass.ss33;
								ss34=nonSphericalClass.ss34;
								qexttxt=nonSphericalClass.cextin;
								qscatxt=nonSphericalClass.csca;
								qabstxt=nonSphericalClass.cabsin;
								albedotxt=nonSphericalClass.walb;
								gtxt=nonSphericalClass.asymm;
								double sum=0;
								for(int i=10;i<171;i++)
								{
									String seq=""+i;
									String mm11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
									d[i]=(Double.valueOf(mm11)-TuScatBean.ss11[i])*(Double.valueOf(mm11)-TuScatBean.ss11[i]);
									sum=sum+d[i];
									System.out.println("hiiiii"+i+"-----"+mm11+"---"+TuScatBean.ss11[i]);
								}
								double rmserror;
								rmserror=Math.sqrt(sum/181);
								if(errorMin==-1.0)
								{
									errorMin=rmserror;
									valueGetrefimasphero=j;
								}
								else if(rmserror<errorMin)
								{
									errorMin=rmserror;
									valueGetrefimasphero=j;
								}
							}
		
							RequestContext context = RequestContext.getCurrentInstance();  
							FacesMessage msg11 = null;
							msg11 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Minimum Error is:", ""+errorMin+"\n"+"Radius"+valueGetrefimasphero);
							FacesContext.getCurrentInstance().addMessage(null, msg11);
							System.out.println("Wave Length Value :"+waveLengthText+"\n");
							System.out.println("Refractive Index(Real) :"+ radiusImgRealTextBox+"\n");
							System.out.println("Refractive Index(Imaginary) :"+ valueGetrefimasphero+"\n");
							System.out.println("Accuracy of T-Matrix Computation :"+accuracyTextBox+"\n");
							System.out.println("A/B :"+a_bOrC_LTextBox+"\n");
							System.out.println("Equivalent Sphere Radius:"+eqvRadiusTextBox+"\n");
		
						}
						/////////////End of Refractive Index and Imaginary part and spheroid
					}  
			/////////////End of Refractive Index and Imaginary part
					
					
					
					
					
					
					
					
					
				}
				///////////////end of Refractive
				
				
				
			
		}
		////end of the try 
		
			
				
				
				
				
			
				
				
				
				
				
		    
			
		
		catch(Exception e)
		{
			e.printStackTrace();	
		}
		
		
      
        
	}
	//////////////////////////////////////For Use with exprimental value page////////////////////////////////////////////
	public void ploatValueExprimentSingle()
	{
		try{
			linearModel = new CartesianChartModel();  
			linearModel2 = new CartesianChartModel();
			linearModel3 = new CartesianChartModel();
			linearModel4 = new CartesianChartModel();
			linearModel1 = new CartesianChartModel();


			  
			        

				LineChartSeries series6 = new LineChartSeries(); 
				LineChartSeries series7 = new LineChartSeries(); 
				LineChartSeries series8 = new LineChartSeries(); 
				LineChartSeries series9 = new LineChartSeries(); 
				LineChartSeries series10 = new LineChartSeries(); 


				try{
						StringReader sr= new StringReader(minimunError);
						BufferedReader br = new BufferedReader(sr);
						ReadWithScanner readWithScanner=new ReadWithScanner();;
						String strLine;
						//Read File Line By Line
						String test="";
						while ((strLine = br.readLine()) != null)   {
						  // Print the content on the console
						  //System.out.println (strLine);
							test+=strLine+"\n";
							readWithScanner.processLine(strLine);
						}
						
				      String ss11 = "";
				      String ss12 = "";
				      String ss33 = "";
				      String ss34 = "";
				
				      
				      for(int i=10;i<171;i++) {
						String seq=""+i;
						ss11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
						ss12=readWithScanner.minimumErrorPojo.getS12().get(seq.trim());
						ss33=readWithScanner.minimumErrorPojo.getS33().get(seq.trim());
						ss34=readWithScanner.minimumErrorPojo.getS34().get(seq.trim());
						series6.set(i,Double.parseDouble(ss11));
						series6.setLabel("ss11");
				        series6.setMarkerStyle("diamond");
				        
						series7.set(i,Double.parseDouble(ss12));
						series7.setLabel("ss12");
				        series7.setMarkerStyle("diamond");
				        
						series8.set(i,Double.parseDouble(ss33));
						series8.setLabel("ss33");
				        series8.setMarkerStyle("diamond");
				        
						series9.set(i,Double.parseDouble(ss34));
						series9.setLabel("ss34");
				        series9.setMarkerStyle("diamond");
				        
						series10.set(i,Math.log10(Double.parseDouble(ss11)));
						series10.setLabel("ss11");
				        series10.setMarkerStyle("diamond");
					}
					
					
					
					
					

						//Close the input stream
						br.close();
						}
						catch (Exception e) {
							// TODO: handle exception
							e.printStackTrace();
						}
			  
			        
			  

					linearModel.addSeries(series6);  
			        linearModel2.addSeries(series7); 
			        linearModel3.addSeries(series8); 
			        linearModel4.addSeries(series9);
			        linearModel1.addSeries(series10);
			       FacesContext context = FacesContext.getCurrentInstance();
			       context = FacesContext.getCurrentInstance();
			       context.getExternalContext().redirect(context.getExternalContext().getRequestContextPath() + "/graph123.jsf");
			}
			catch(Exception e)
			{
			e.printStackTrace();
			}
			        
			       // return "graph.xhtml?faces-redirect=true";
	}
	//////////////////////////////For Use with theory page only/////////////////////////////////
	public void ploatValueExpriment()
	{

	try{
	linearModel = new CartesianChartModel();  
	linearModel2 = new CartesianChartModel();
	linearModel3 = new CartesianChartModel();
	linearModel4 = new CartesianChartModel();
	linearModel1 = new CartesianChartModel();


	  
	        LineChartSeries series1 = new LineChartSeries();  
	        series1.setLabel("ss11(Linear)");
	        series1.setMarkerStyle("diamond");
	        for(int i=10;i<171;i++)
	        {
	        series1.set(i,ss11[i]);  
	        }
	  
	        LineChartSeries series5 = new LineChartSeries();  
	        series5.setLabel("ss11(Log)"); 
	        series5.setMarkerStyle("diamond");
	        for(int i=10;i<171;i++)
	        {
	        series5.set(i,Math.log10(ss11[i]));  
	        }
	       
	  
	        LineChartSeries series2 = new LineChartSeries();  
	        series2.setLabel("ss12");  
	        series2.setMarkerStyle("diamond");  
	        
	        for(int i=10;i<171;i++)
	        {
	        series2.set(i,ss12[i]);  
	        }
	        
	        LineChartSeries series3 = new LineChartSeries();  
	        series3.setLabel("ss33");
	        series3.setMarkerStyle("diamond");
	        for(int i=10;i<171;i++)
	        {
	        series3.set(i,ss33[i]);  
	        }
	        LineChartSeries series4 = new LineChartSeries();  
	        series4.setLabel("ss34");  
	        series4.setMarkerStyle("diamond");  
	        
	        for(int i=10;i<171;i++)
	        {
	        series4.set(i,ss34[i]);  
	        }

		LineChartSeries series6 = new LineChartSeries(); 
		LineChartSeries series7 = new LineChartSeries(); 
		LineChartSeries series8 = new LineChartSeries(); 
		LineChartSeries series9 = new LineChartSeries(); 
		LineChartSeries series10 = new LineChartSeries(); 


		try{
				StringReader sr= new StringReader(minimunError);
				BufferedReader br = new BufferedReader(sr);
				ReadWithScanner readWithScanner=new ReadWithScanner();;
				String strLine;
				//Read File Line By Line
				String test="";
				while ((strLine = br.readLine()) != null)   {
				  // Print the content on the console
				  //System.out.println (strLine);
					test+=strLine+"\n";
					readWithScanner.processLine(strLine);
					
					//System.out.println(""+strLine);
				}
				
				
		      String ss11 = "";
		      String ss12 = "";
		      String ss33 = "";
		      String ss34 = "";
		
		      
		      for(int i=10;i<171;i++) {
				String seq=""+i;
				ss11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
				ss12=readWithScanner.minimumErrorPojo.getS12().get(seq.trim());
				ss33=readWithScanner.minimumErrorPojo.getS33().get(seq.trim());
				ss34=readWithScanner.minimumErrorPojo.getS34().get(seq.trim());
				series6.set(i,Double.parseDouble(ss11));
				series6.setLabel("ss11");  
		        series6.setMarkerStyle("diamond");
		        
				series7.set(i,Double.parseDouble(ss12));
				series7.setLabel("ss12");  
		        series7.setMarkerStyle("diamond");
		        
				series8.set(i,Double.parseDouble(ss33));
				series8.setLabel("ss33");  
		        series8.setMarkerStyle("diamond");
		        
				series9.set(i,Double.parseDouble(ss34));
				series9.setLabel("ss34");  
		        series9.setMarkerStyle("diamond");
		        
				series10.set(i,Math.log10(Double.parseDouble(ss11)));
				series10.setLabel("ss11");  
		        series10.setMarkerStyle("diamond");
				
			}
			/*int sum=0;
			for(double j=minRadius;j<=maxRadius;j=j+stepRadius)
			{
		main_Function();
			for (int i=0;i<181;i++) {
				String seq=""+i;
				
				ss11=readWithScanner.minimumErrorPojo.getS11().get(seq.trim());
				
				sum+=((Double.valueOf(ss11)-(TuScatBean.ss11[i]))*(ss11Com[i]-(TuScatBean.ss11[i])));
			}
			double rmserror;
			rmserror=Math.sqrt(sum/181);
			if(errorMin==-1.0)
			{
			errorMin=rmserror;
			valueGet=j;
			}
			else if(rmserror<errorMin)
			{
			errorMin=rmserror;
			valueGet=j;
			}
			}*/
				//Close the input stream
				br.close();
				}
				catch (Exception e) {
					// TODO: handle exception
					e.printStackTrace();
				}
	  
	        
	  
	        linearModel.addSeries(series1);  
	        linearModel2.addSeries(series2); 
	        linearModel3.addSeries(series3); 
	        linearModel4.addSeries(series4);
	        linearModel1.addSeries(series5);
	        linearModel.addSeries(series6);  
	        linearModel2.addSeries(series7); 
	        linearModel3.addSeries(series8); 
	        linearModel4.addSeries(series9);
	        linearModel1.addSeries(series10);
	       FacesContext context = FacesContext.getCurrentInstance();
	       context = FacesContext.getCurrentInstance();
	       context = FacesContext.getCurrentInstance();
	       context.getExternalContext().redirect(context.getExternalContext().getRequestContextPath() + "/graph123.jsf");
	}
	catch(Exception e)
	{
	e.printStackTrace();
	}
	        
	       // return "graph.xhtml?faces-redirect=true";
	}
	
	
	
	public void reset()
	{
		
		
		
		
		 sizeRadioButton=false;
		
		 refractiveRadioButton=false;
		
		
		 realindexRadioButton=false;
		 imagindexRadioButton=false;
	
		 sphereRadioButton=false;
		 spheroidRadioButton=false;
		 cyclinderRadioButton=false;
		
		
		
		 boolsizeRI=false;
		 boolsizeRISpCy=false;
		 boolsizeRISproid=false;
		 boolsizeABCD=false;
		 
		 boolrefractive=false;
		 boolrefPara=false;
		
		
		
		
		
		 fixedrrefractiveIndexTextBox="";
		 fixedirefractiveIndexTextBox="";
		
		
		
		
		 selectsizeRef="";
		
		
		 minRadiusTextBox="";
		 maxRadiusTextBox="";
		 stepRadiusTextBox="";
		
		
		 rrefractiveIndexTextBox="";
		 irefractiveIndexTextBox="";
		 waveLengthText="";
		 selectcomboTypesShap="";
		 combosizeDis="";
		 selectcombosizeDis="";
		 radiusTextBox="";
		 lowestPartRadiusText="";
		 highestPartRadiusText="";
		 sigmaText="";
		 modalradiusText="";
		 extinctionCoefficentValue="";
		 scattereingCoefficentValue="";
		 absorptionCoefficentValue="";
		 singleScatteringValue="";
		 assymmetricLabeltValue="";
		
		
		
		minRadiusTextBox="";
		maxRadiusTextBox="";
		stepRadiusTextBox="";
		accuracyTextBox="";
		// a_bOrC_LTextBox;
		 selectcomboTypesPara="";
		 eqvRadiusTextBox="";
		 minAccuracyTextBox="";
		 maxAccuracyTextBox="";
		 stepAccuracyTextBox="";
		 displayTextABVE="";
		 displayTextComBo="";
		 displayTextRealImag="";
		 displayTextLegen="";
		 displayTextLast="";
		
		 displayexpectedRange="";
		
		 selectcomboTypesShapImg="";
		
		
		  radiusImgRealTextBox="";
		  minirefractiveIndexTextBox="";
		  maxrefractiveIndexTextBox="";
		  steprefractiveIndexTextBox="";
		  a_bOrC_LTextBox="";
	}
	
	
	
	public void onChangeGeometry()
	{
		
		//System.out.println("test test");
		if(selectcomboTypesShap.equalsIgnoreCase("Select"))
		{
			boolsizeRI=false;
			boolsizeRISpCy=false;
			boolsizeRISproid=false;
			boolsizeABCD1=false;
			boolsizeABCD=false;
			boolsizeABCDspheroidab=false;
			boolsizeABCDspheroiddl=false;
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Sphere"))
		{
			boolsizeRI=true;
			boolsizeRISpCy=false;
			boolsizeRISproid=false;
			boolsizeABCD1=false;
			boolsizeABCD=false;
			boolsizeABCDspheroidab=false;
			boolsizeABCDspheroiddl=false;
			
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Cylinder"))
		{
			boolsizeRI=false;	
			boolsizeRISpCy=true;
			boolsizeRISproid=false;
			boolsizeABCDspheroidab=false;
			boolsizeABCDspheroiddl=false;
			displayTextComBo="D/L Ratio";
			displayTextComBoEqui="Equivalent Radius";
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Spheroid"))
		{
			boolsizeRI=false;
			boolsizeRISpCy=false;
			boolsizeRISproid=true;
			displayTextComBo="A/B Ratio";
			displayTextComBoEqui="Equivalent Radius";
			boolsizeABCD1=false;
			boolsizeABCD=false;
		}
	}
	public void onChangeTypesParameter()
	{
		if(selectParameter.equalsIgnoreCase("Select"))
		{
			boolsizeABCD1=false;
			boolsizeABCD=false;
			boolsizeABCDspheroidab=false;
			boolsizeABCDspheroiddl=false;
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Cylinder")&&selectParameter.equalsIgnoreCase("D/L Ratio"))
		{
			boolsizeABCD=true;
			boolsizeABCD1=false;
			boolsizeABCDspheroidab=false;
			boolsizeABCDspheroiddl=false;
			displayTextABVE="Volume Equivalent Radius";
			displayexpectedRange="Expected Range of  D/L Ratio";
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Cylinder")&&selectParameter.equalsIgnoreCase("Equivalent Radius"))
		{
			displayTextABVE="D/L Ratio";
			displayexpectedRange="Expected Range of  Equivalent Radius";
			boolsizeABCD1=true;
			boolsizeABCD=false;
			boolsizeABCDspheroidab=false;
			boolsizeABCDspheroiddl=false;
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Spheroid")&&selectParameter.equalsIgnoreCase("A/B Ratio"))
		{
			displayTextABVE="Volume Equivalent Radius";
			displayexpectedRange="Expected Range of  A/B Ratio";
			boolsizeABCDspheroidab=true;
			boolsizeABCDspheroiddl=false;
			boolsizeABCD1=false;
			boolsizeABCD=false;
		}
		else if(selectcomboTypesShap.equalsIgnoreCase("Spheroid")&&selectParameter.equalsIgnoreCase("Equivalent Radius"))
		{
			displayTextABVE="A/B Ratio";
			displayexpectedRange="Expected Range of  Equivalent Radius";
			boolsizeABCDspheroiddl=true;
			boolsizeABCDspheroidab=false;
			boolsizeABCD1=false;
			boolsizeABCD=false;
		}
		
	}
	
	
	
	
	
	
	
	  
	  
	  
	 public CartesianChartModel getLinearModel() {  
	        return linearModel;  
	    }  
	 
	 
	
	
	public void checkSize()
	{
		//sizeRadioButton=false;
		//refractiveRadioButton=false;
		boolsizeABCD=false;
		boolsizeRISpCy=false;
		boolsizeRISproid=false;
		boolrefrctiveSphere=false;
		boolrefrctiveCylinder=false;
		boolrefrctiveSpheroid=false;
		boolsizeRI=false;
		boolrefPara=false;
		boolrefractive=false;
		boolrefractive=false;
		
		if(sizeRadioButton)
		{
			sizeRadioButton=true;
			refractiveRadioButton=false;
		}
		else if(refractiveRadioButton)
		{
			sizeRadioButton=false;
			refractiveRadioButton=true;
		}
	}  
		
	
	
	public void checkRealImaginary()
	{
		//imagindexRadioButton=false;
		//realindexRadioButton=false;
		boolrefrctiveSphere=false;
		boolrefrctiveCylinder=false;
		boolrefrctiveSpheroid=false;
		boolsizeRI=false;
		boolrefPara=false;
		boolrefractive=false;
		boolrefractive=false;
		if(realindexRadioButton)
		{
			realindexRadioButton=true;
			imagindexRadioButton=false;
			refractiveRadioButton=true;
			boolrefractive=true;
			displayTextRealImag="Enter the Imaginary part :";
			displayTextLegen="Expected Range of Real Part of the Refractive Index"+"\n";
		}
		else if(imagindexRadioButton)
		{
			imagindexRadioButton=true;
			realindexRadioButton=false;
			refractiveRadioButton=true;
			boolrefractive=true;
			displayTextRealImag="Enter the real part :";
			displayTextLegen="Expected Range of Imaginary Part of the Refractive Index";
		}
	
	}
	
	public void onChangeGeometryImg()
	{
		//System.out.println("test test");
		if(selectcomboTypesShapImg.equalsIgnoreCase("Select"))
		{
			boolrefrctiveSphere=false;
			boolrefrctiveCylinder=false;
			boolrefrctiveSpheroid=false;
		}
		else if(selectcomboTypesShapImg.equalsIgnoreCase("Sphere"))
		{
			boolrefrctiveSphere=true;
			boolrefrctiveCylinder=false;
			boolrefrctiveSpheroid=false;
		}
		else if(selectcomboTypesShapImg.equalsIgnoreCase("Cylinder"))
		{
			boolrefrctiveSphere=false;
			boolrefrctiveCylinder=true;
			boolrefrctiveSpheroid=false;
			displayTextLast="Diameter to Length Ratio,D/L";
		}
		else if(selectcomboTypesShapImg.equalsIgnoreCase("Spheroid"))
		{
			boolrefrctiveSphere=false;
			boolrefrctiveCylinder=false;
			boolrefrctiveSpheroid=true;
			displayTextLast="Horizontal to Rotational Axis Ratio,A/B";
		}
		
	}
	/*public void  checkGeo()
	{
		if(sphereRadioButton)
		{
			cyclinderRadioButton=false;
			spheroidRadioButton=false;
		}
		if(spheroidRadioButton)
		{
			sphereRadioButton=false;
			cyclinderRadioButton=false;
		}
		if(cyclinderRadioButton)
		{
			sphereRadioButton=false;
			spheroidRadioButton=false;
			
		}
	}*/
	
	public void onChangecombosizeDis()
	{
		//System.out.println("selectcomboTypesShap"+selectcomboTypesShap);
		if(null!=selectcombosizeDis&&!selectcombosizeDis.equalsIgnoreCase("Select"))
		{
			if(selectcombosizeDis.equalsIgnoreCase("Monodisperse"))
			{
				boolshape=false;
				boolshapenot=true;
				//System.out.println("selectcomboTypesShap"+selectcomboTypesShap);
			}
			else
			{
				boolshape=true;
				boolshapenot=false;
				//System.out.println("selectcomboTypesShap"+selectcomboTypesShap);
			}
			if(selectcombosizeDis.equalsIgnoreCase("Gamma"))
			{
				labelalfasigma="Alfa";
			}
			else
			{
				labelalfasigma="Sigma";
			}
			
		}
	}
	
	
	public void onChangecomboTypesShap()
	{
		if(null!=selectcomboTypesShap&&!selectcomboTypesShap.equalsIgnoreCase("Select"))
		{
			if(selectcomboTypesShap.equalsIgnoreCase("Sphere"))
			{
				
				boolsphere=false;
				//System.out.println("selectcomboTypesShap"+selectcomboTypesShap);
			}
			else
			{
				boolsphere=true;
				
				if(selectcomboTypesShap.equalsIgnoreCase("Cylinder"))
				{
				np=-2;
				}
				else
				{
				np=-1;
				}
				
				//System.out.println("selectcomboTypesShap"+selectcomboTypesShap);
			}
			
			
		}
		
	}
	
	
	public void ploatValue()
	{
		
		try{
		linearModel = new CartesianChartModel();  
		linearModel2 = new CartesianChartModel();
		linearModel3 = new CartesianChartModel();
		linearModel4 = new CartesianChartModel();
		linearModel1 = new CartesianChartModel();
		
		
		  
        LineChartSeries series1 = new LineChartSeries();  
        series1.setLabel("ss11(Linear)");  
        for(int i=0;i<181;i++)
        {
        	 series1.set(i,ss11[i]);  
        }
  
        LineChartSeries series5 = new LineChartSeries();  
        series5.setLabel("ss11(Log)");  
        for(int i=0;i<181;i++)
        {
        	 series5.set(i,Math.log10(ss11[i]));  
        }
       
  
        LineChartSeries series2 = new LineChartSeries();  
        series2.setLabel("ss12");  
        series2.setMarkerStyle("diamond");  
        
        for(int i=0;i<181;i++)
        {
        	 series2.set(i,ss12[i]);  
        }
        
        LineChartSeries series3 = new LineChartSeries();  
        series3.setLabel("ss33");  
        for(int i=0;i<181;i++)
        {
        	 series3.set(i,ss33[i]);  
        }
        LineChartSeries series4 = new LineChartSeries();  
        series4.setLabel("ss34");  
        series4.setMarkerStyle("diamond");  
        
        for(int i=0;i<181;i++)
        {
        	 series4.set(i,ss34[i]);  
        }
  
        
  
        linearModel.addSeries(series1);  
        linearModel2.addSeries(series2); 
        linearModel3.addSeries(series3); 
        linearModel4.addSeries(series4);
        linearModel1.addSeries(series5);
       FacesContext context = FacesContext.getCurrentInstance();
       context = FacesContext.getCurrentInstance();
       context = FacesContext.getCurrentInstance();
       context.getExternalContext().redirect(context.getExternalContext().getRequestContextPath() + "/graph123.jsf");
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
        
       // return "graph.xhtml?faces-redirect=true";
	}
	public void saveTheorialDataValue()
	{

		try
		{
			int m=0;
			String filename="../webapps/WEBSCATwebapplication/New Record/TUSCAT.xls";
			HSSFWorkbook hwb=new HSSFWorkbook();
			HSSFSheet sheet =  hwb.createSheet("new sheet");
			
			
			for(int i=0;i<ss11.length;i++)
			{
			DecimalFormat number = new DecimalFormat("#.000000");
			HSSFRow row=sheet.createRow((short)i);
			row.createCell((short) 0).setCellValue(""+i);
			row.createCell((short) 1).setCellValue(","+number.format(ss11[i]));
			row.createCell((short) 2).setCellValue(","+number.format(ss12[i]));
			row.createCell((short) 3).setCellValue(","+number.format(ss33[i]));
			row.createCell((short) 4).setCellValue(","+number.format(ss34[i]));
			
			
			//System.out.println(""+number.format(ss11[i]));
			//System.out.println(".."+ss11[i]+".."+ss12[i]+".."+ss33[i]+".."+ss34[i]);
			
			
			}
			
			
			
			
			
			FileOutputStream fileOut =  new FileOutputStream(filename);
			
			hwb.write(fileOut);
			fileOut.close();
			InputStream stream = ((ServletContext)FacesContext.getCurrentInstance().getExternalContext().getContext()).getResourceAsStream("/New Record/TUSCAT.xls");  
	        file = new DefaultStreamedContent(stream, "New Record/xls", "TUSCAT.xls");
			
			System.out.println("Your File Has been written!");
		}	
	      catch(Exception e)
	      {
	    	  e.printStackTrace();
	      }
	        
	    
	}
	
	
	
	public String getLowestPartRadiusText() {
		return lowestPartRadiusText;
	}
	public void setLowestPartRadiusText(String lowestPartRadiusText) {
		this.lowestPartRadiusText = lowestPartRadiusText;
	}
	public String getHighestPartRadiusText() {
		return highestPartRadiusText;
	}
	public void setHighestPartRadiusText(String highestPartRadiusText) {
		this.highestPartRadiusText = highestPartRadiusText;
	}
	public String getSigmaText() {
		return sigmaText;
	}
	public void setSigmaText(String sigmaText) {
		this.sigmaText = sigmaText;
	}
	public String getModalradiusText() {
		return modalradiusText;
	}
	public void setModalradiusText(String modalradiusText) {
		this.modalradiusText = modalradiusText;
	}
	public String getSelectcomboTypesShap() {
		return selectcomboTypesShap;
	}
	public void setSelectcomboTypesShap(String selectcomboTypesShap) {
		this.selectcomboTypesShap = selectcomboTypesShap;
	}
	public String getRrefractiveIndexTextBox() {
		return rrefractiveIndexTextBox;
	}
	public void setRrefractiveIndexTextBox(String rrefractiveIndexTextBox) {
		this.rrefractiveIndexTextBox = rrefractiveIndexTextBox;
	}
	public String getIrefractiveIndexTextBox() {
		return irefractiveIndexTextBox;
	}
	public void setIrefractiveIndexTextBox(String irefractiveIndexTextBox) {
		this.irefractiveIndexTextBox = irefractiveIndexTextBox;
	}
	public String getWaveLengthText() {
		return waveLengthText;
	}
	public void setWaveLengthText(String waveLengthText) {
		this.waveLengthText = waveLengthText;
	}
	public String getCombosizeDis() {
		return combosizeDis;
	}
	public void setCombosizeDis(String combosizeDis) {
		this.combosizeDis = combosizeDis;
	}
	public String getSelectcombosizeDis() {
		return selectcombosizeDis;
	}
	public void setSelectcombosizeDis(String selectcombosizeDis) {
		this.selectcombosizeDis = selectcombosizeDis;
	}
	public String getAccuracyTextBox() {
		return accuracyTextBox;
	}
	public void setAccuracyTextBox(String accuracyTextBox) {
		this.accuracyTextBox = accuracyTextBox;
	}
	public String getA_bOrC_LTextBox() {
		return a_bOrC_LTextBox;
	}
	public void setA_bOrC_LTextBox(String orC_LTextBox) {
		a_bOrC_LTextBox = orC_LTextBox;
	}
	public String getExtinctionCoefficentValue() {
		return extinctionCoefficentValue;
	}
	public void setExtinctionCoefficentValue(String extinctionCoefficentValue) {
		this.extinctionCoefficentValue = extinctionCoefficentValue;
	}
	public String getScattereingCoefficentValue() {
		return scattereingCoefficentValue;
	}
	public void setScattereingCoefficentValue(String scattereingCoefficentValue) {
		this.scattereingCoefficentValue = scattereingCoefficentValue;
	}
	public String getAbsorptionCoefficentValue() {
		return absorptionCoefficentValue;
	}
	public void setAbsorptionCoefficentValue(String absorptionCoefficentValue) {
		this.absorptionCoefficentValue = absorptionCoefficentValue;
	}
	public String getSingleScatteringValue() {
		return singleScatteringValue;
	}
	public void setSingleScatteringValue(String singleScatteringValue) {
		this.singleScatteringValue = singleScatteringValue;
	}
	public String getAssymmetricLabeltValue() {
		return assymmetricLabeltValue;
	}
	public void setAssymmetricLabeltValue(String assymmetricLabeltValue) {
		this.assymmetricLabeltValue = assymmetricLabeltValue;
	}
	public boolean isBoolsphere() {
		return boolsphere;
	}
	public void setBoolsphere(boolean boolsphere) {
		this.boolsphere = boolsphere;
	}



	public String getRadiusTextBox() {
		return radiusTextBox;
	}



	public void setRadiusTextBox(String radiusTextBox) {
		this.radiusTextBox = radiusTextBox;
	}


	public boolean isBoolshape() {
		return boolshape;
	}


	public void setBoolshape(boolean boolshape) {
		this.boolshape = boolshape;
	}


	public String getLabelalfasigma() {
		return labelalfasigma;
	}


	public void setLabelalfasigma(String labelalfasigma) {
		this.labelalfasigma = labelalfasigma;
	}


	public boolean isBoolshapenot() {
		return boolshapenot;
	}


	public void setBoolshapenot(boolean boolshapenot) {
		this.boolshapenot = boolshapenot;
	}


	public boolean isSizeRadioButton() {
		return sizeRadioButton;
	}


	public void setSizeRadioButton(boolean sizeRadioButton) {
		this.sizeRadioButton = sizeRadioButton;
	}


	public boolean isRefractiveRadioButton() {
		return refractiveRadioButton;
	}


	public void setRefractiveRadioButton(boolean refractiveRadioButton) {
		this.refractiveRadioButton = refractiveRadioButton;
	}
	
	public void calculateValue()
	{
		 System.out.println("Start");
		 //System.out.println("irefractiveIndexTextBox"+irefractiveIndexTextBox);
		 //System.out.println("rrefractiveIndexTextBox"+rrefractiveIndexTextBox);
		 //System.out.println("waveLengthText"+waveLengthText);
		 //System.out.println("lowestPartRadiusText"+lowestPartRadiusText);
		 //System.out.println("highestPartRadiusText"+highestPartRadiusText);
		 //System.out.println("sigmaText"+sigmaText);
		 //System.out.println("modalradiusText"+modalradiusText);
		
		 
		 
		 
		double radlk=0.0;
		double radhk=0.0;
		double stepk;
		int ntotk;
		
		for(int i=0;i<ss11.length;i++)
		{
		ss11[i]=0.0;
		ss12[i]=0.0;
		ss33[i]=0.0;
		ss34[i]=0.0;
		
		}
		if(null!=rrefractiveIndexTextBox&&(rrefractiveIndexTextBox.trim().equalsIgnoreCase("")))
		{
			//FacesContext.getCurrentInstance().addMessage(null, new FacesMessage(FacesMessage.SEVERITY_ERROR,"Sample error message", "PrimeFaces makes no mistakes")); 
		}
		if(null!=irefractiveIndexTextBox&&(irefractiveIndexTextBox.trim().equalsIgnoreCase("")))
		{
			//FacesContext.getCurrentInstance().addMessage(null, new FacesMessage(FacesMessage.SEVERITY_ERROR,"Sample error message", "PrimeFaces makes no mistakes")); 
		}
		if(null!=waveLengthText&&(waveLengthText.trim().equalsIgnoreCase("")))
		{
			//FacesContext.getCurrentInstance().addMessage(null, new FacesMessage(FacesMessage.SEVERITY_ERROR,"Sample error message", "PrimeFaces makes no mistakes")); 
		}
		String strrefre=rrefractiveIndexTextBox;
		refre=Double.parseDouble(strrefre);
		String strrefim=irefractiveIndexTextBox;
		refim=Double.parseDouble(strrefim);
		String strrefmedreal="1";
		double real=Double.parseDouble(strrefmedreal);
		refmed=new Complex(real,0.0);
		refrelk=new Complex(refre,refim);
		ntotk=1; 
		stepk=0.05;
		String strxk1=waveLengthText;
		double lam=Double.parseDouble(strxk1);
		
		if(selectcombosizeDis.equalsIgnoreCase("Monodisperse"))
		{
		if(null!=radiusTextBox&&(radiusTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value forRadius");
		return;
		}
		
		String strxk2 =radiusTextBox;
		xk=(2.0*PI*Double.parseDouble(strxk2)*real)/Double.parseDouble(strxk1);
		//System.out.println("xk"+xk);

		if(selectcomboTypesShap.equalsIgnoreCase("Sphere"))
		{

		/**********************************/
		main_Function();
		}
		else
		{
		if(null!=a_bOrC_LTextBox&&a_bOrC_LTextBox.trim().equalsIgnoreCase(""))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for A/B Or C/L Field");
		return;
		}
		if(null!=accuracyTextBox&&(accuracyTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Accuracy Field");
		return;
		}
		double axi=Double.parseDouble(strxk2);
		String streps=a_bOrC_LTextBox;
		double eps=Double.parseDouble(streps);
		String strddelt=accuracyTextBox;
		double ddelt=Double.parseDouble(strddelt);
		radlk=0.9999999*axi;
		radhk=1.0000001*axi;
		String ret=nonSphericalClass.nonsp(refre,refim,1,axi,radlk,radhk,7,eps,ddelt,-1,np,lam);
		//System.out.println("");
		if(null!=ret&&!ret.equalsIgnoreCase(""))
		{
		//JOptionPane.showMessageDialog(null, ret);
		return;
		}
		ss11=nonSphericalClass.ss11;
		ss12=nonSphericalClass.ss12;
		ss33=nonSphericalClass.ss33;
		ss34=nonSphericalClass.ss34;
		qexttxt=nonSphericalClass.cextin;
		qscatxt=nonSphericalClass.cscat;
		qabstxt=nonSphericalClass.cabsin;
		albedotxt=nonSphericalClass.walb;
		gtxt=nonSphericalClass.asymm;
		
		}
		}
		else if(selectcombosizeDis.equalsIgnoreCase("Gamma"))
		{
		if(null!=lowestPartRadiusText&&(null==lowestPartRadiusText ||lowestPartRadiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Highest Partical Radius Field");
		return;
		}
		if(null!=highestPartRadiusText&&(null==highestPartRadiusText||highestPartRadiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Lowest Partical Radius Field");
		return;
		}
		if(null!=modalradiusText&&(null==modalradiusText||modalradiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Modal Radius Field");
		return;
		}
		if(null!=sigmaText&&(null==sigmaText||sigmaText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Sigma Field");
		return;
		}
		String strradlk=lowestPartRadiusText;
		radlk=Double.parseDouble(strradlk);
		String strradhk=highestPartRadiusText;
		radhk=Double.parseDouble(strradhk);
		wavel=Double.parseDouble(strxk1);
		String strxk2=modalradiusText;
		String stralfa=sigmaText;
		double rck=Double.parseDouble(strxk2);
		double alfak=Double.parseDouble(stralfa);

		//System.out.println("//"+radlk+"//"+radhk+"//"+wavel+"//"+rck+"//"+alfak);
		if(selectcomboTypesShap.equalsIgnoreCase("Sphere"))
		{
		extra(stepk,radlk,radhk,ntotk,rck,alfak,3);
		}
		else
		{
		if(null!=a_bOrC_LTextBox&&(null==a_bOrC_LTextBox||a_bOrC_LTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for A/B Or C/L Field");
		return;
		}
		if(null!=accuracyTextBox&&(null==accuracyTextBox||accuracyTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Accuracy Field");
		return;
		}
		String streps=a_bOrC_LTextBox;
		double eps=Double.parseDouble(streps);
		String strddelt=accuracyTextBox;
		double ddelt=Double.parseDouble(strddelt);
		
		
		String ret=nonSphericalClass.nonsp(refre,refim,1,rck,radlk,radhk,alfak,eps,ddelt,5,np,wavel);
		
		//System.out.println("////"+radlk+"////"+radhk+"////"+wavel+"////"+rck+"////"+alfak);
		if(null!=ret&&!ret.equalsIgnoreCase(""))
		{
		//JOptionPane.showMessageDialog(null, ret);
		return;
		}
		ss11=nonSphericalClass.ss11;
		ss12=nonSphericalClass.ss12;
		ss33=nonSphericalClass.ss33;
		ss34=nonSphericalClass.ss34;
		qexttxt=nonSphericalClass.cextin;
		qscatxt=nonSphericalClass.cscat;
		qabstxt=nonSphericalClass.cabsin;
		albedotxt=nonSphericalClass.walb;
		gtxt=nonSphericalClass.asymm;
		}
		}
		else if(selectcombosizeDis.equalsIgnoreCase("Normal"))
		{
		if(null!=lowestPartRadiusText&&(null==lowestPartRadiusText||lowestPartRadiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Highest Partical Radius Field");
		return;
		}
		if(null!=highestPartRadiusText&&(null==highestPartRadiusText||highestPartRadiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Lowest Partical Radius Field");
		return;
		}
		if(null!=modalradiusText&&(null==modalradiusText||modalradiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Modal Radius Field");
		return;
		}
		if(null!=sigmaText&&(null==sigmaText||sigmaText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Sigma Field");
		return;
		}
		String strradlk=lowestPartRadiusText;
		radlk=Double.parseDouble(strradlk);
		String strradhk=highestPartRadiusText;
		radhk=Double.parseDouble(strradhk);
		 
		wavel=Double.parseDouble(strxk1);
		String strxk2=modalradiusText;
		String stralfa=sigmaText;
		double rck=Double.parseDouble(strxk2);
		double alfak=Double.parseDouble(stralfa);
		if(selectcomboTypesShap.equalsIgnoreCase("Sphere"))
		{
		extra(stepk,radlk,radhk,ntotk,rck,alfak,3);
		}
		else
		{
		if(null!=a_bOrC_LTextBox&&(null==a_bOrC_LTextBox||a_bOrC_LTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for A/B Or C/L Field");
		return;
		}
		if(null!=accuracyTextBox&&(null==accuracyTextBox||accuracyTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Accuracy Field");
		return;
		}
		String streps=a_bOrC_LTextBox;
		double eps=Double.parseDouble(streps);
		String strddelt=accuracyTextBox;
		double ddelt=Double.parseDouble(strddelt);
		String ret=nonSphericalClass.nonsp(refre,refim,2,rck,radlk,radhk,alfak,eps,ddelt,5,np,wavel);
		if(null!=ret&&!ret.equalsIgnoreCase(""))
		{
		//JOptionPane.showMessageDialog(null, ret);
		return;
		}
		ss11=nonSphericalClass.ss11;
		ss12=nonSphericalClass.ss12;
		ss33=nonSphericalClass.ss33;
		ss34=nonSphericalClass.ss34;
		qexttxt=nonSphericalClass.cextin;
		qscatxt=nonSphericalClass.cscat;
		qabstxt=nonSphericalClass.cabsin;
		albedotxt=nonSphericalClass.walb;
		gtxt=nonSphericalClass.asymm;
		}
		}
		else if(selectcombosizeDis.equalsIgnoreCase("Lognormal"))
		{
		if(null!=lowestPartRadiusText&&(null==lowestPartRadiusText||lowestPartRadiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Highest Partical Radius Field");
		return;
		}
		if(null!=highestPartRadiusText&&(null==highestPartRadiusText||highestPartRadiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Lowest Partical Radius Field");
		return;
		}
		if(null!=modalradiusText&&(null==modalradiusText||modalradiusText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Modal Radius Field");
		return;
		}
		if(null!=sigmaText&&(null==sigmaText||sigmaText.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Sigma Field");
		return;
		}
		String strradlk=lowestPartRadiusText;
		radlk=Double.parseDouble(strradlk);
		String strradhk=highestPartRadiusText;
		radhk=Double.parseDouble(strradhk);
		wavel=Double.parseDouble(strxk1);
		String strxk2=modalradiusText;
		String stralfa=sigmaText;
		double rck=Double.parseDouble(strxk2);
		double alfak=Double.parseDouble(stralfa);
		if(selectcomboTypesShap.equalsIgnoreCase("Sphere"))
		{
		extra(stepk,radlk,radhk,ntotk,rck,alfak,3);
		}
		else
		{
		if(null!=a_bOrC_LTextBox&&(null==a_bOrC_LTextBox||a_bOrC_LTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for A/B Or C/L Field");
		return;
		}
		if(null!=accuracyTextBox&&(null==accuracyTextBox||accuracyTextBox.trim().equalsIgnoreCase("")))
		{
		//JOptionPane.showMessageDialog(null, "Please Enter The Value for Accuracy Field");
		return;
		}
		String streps=a_bOrC_LTextBox;
		double eps=Double.parseDouble(streps);
		String strddelt=accuracyTextBox;
		double ddelt=Double.parseDouble(strddelt);
		String ret=nonSphericalClass.nonsp(refre,refim,3,rck,radlk,radhk,alfak,eps,ddelt,5,np,wavel);
		if(null!=ret&&!ret.equalsIgnoreCase(""))
		{
		//JOptionPane.showMessageDialog(null, ret);
		return;
		}
		ss11=nonSphericalClass.ss11;
		ss12=nonSphericalClass.ss12;
		ss33=nonSphericalClass.ss33;
		ss34=nonSphericalClass.ss34;
		qexttxt=nonSphericalClass.cextin;
		qscatxt=nonSphericalClass.cscat;
		qabstxt=nonSphericalClass.cabsin;
		albedotxt=nonSphericalClass.walb;
		gtxt=nonSphericalClass.asymm;
		}
		if(qexttxt!=0.0)
		{
		//boolCalculation=true;
		}
		
			
		}
		listofdata.clear();
		DecimalFormat df = new DecimalFormat("#0.0000");
		for(int i=0;i<ss11.length;i++)
		{
			listofdata.add(new TuScatPojo(i, ss11[i], ss12[i], ss33[i], ss34[i]));
		}
		
		
		extinctionCoefficentValue=df.format(qexttxt);
		scattereingCoefficentValue=df.format(qscatxt);
		absorptionCoefficentValue=df.format(qabstxt);
		singleScatteringValue=df.format(albedotxt);
		assymmetricLabeltValue=df.format(gtxt);
		
		//System.out.println("extinctionCoefficentValue"+extinctionCoefficentValue);
		//System.out.println("scattereingCoefficentValue"+scattereingCoefficentValue);
		//System.out.println("absorptionCoefficentValue"+absorptionCoefficentValue);
		//System.out.println("singleScatteringValue"+singleScatteringValue);
		//System.out.println("assymmetricLabeltValue"+assymmetricLabeltValue);
		System.out.println("End");
		
		
		

	}
	 public static void main_Function()
	 {
	 int count=0;
	 refrelk=new Complex(refre,refim);
	 refrelk=refrelk.div(refmed);
	 dang=(PI/2.0)/(double)(nangk-1);
	 bhmie();
	 gtxt=gscatxt/qscatxt;
	 qprtxt=qexttxt-gscatxt;
	 qabstxt=qexttxt-qscatxt;
	 albedotxt=qscatxt/qexttxt;
	 //System.out.println("gtxt"+gtxt);
	 //System.out.println("qprtxt"+qprtxt);
	 //System.out.println("qabstxt"+qabstxt);
	 //System.out.println("albedotxt"+albedotxt);
	 s11nor=(0.5)*(Math.pow(s2k[1].mod(),2.0)+Math.pow(s1k[1].mod(),2.0));
	 nan=2*nangk-1;
	 for(j=1;j<=nan;j++)
	 {
	 s11=(0.5)*(Math.pow(s2k[j].mod(),2.0)+Math.pow(s1k[j].mod(),2.0));
	 s12=(0.5)*(Math.pow(s2k[j].mod(),2.0)-Math.pow(s1k[j].mod(),2.0));
	// System.out.println("s11"+s11);
	// System.out.println("s1k[j]"+s1k[j]);
	// System.out.println("s2k[j]"+s2k[j]);
	 pol=-s12/s11;
	 s33=(s2k[j].times(s1k[j].conj())).real();
	 s33=s33/s11;
	 s34=(s2k[j].times(s1k[j].conj())).imag();
	 s34=s34/s11;
	 s11=s11/s11nor;
	 ang=dang*(j-1.0)*(180.0/PI);
	 ss11[count]=s11;
	 ss12[count]=pol;
	 ss33[count]=s33;
	 ss34[count]=s34;
	 count++;
	 }
	 }

	 /***********************************************************************/
	 /******************BHMIE function for Mie calculations******************/
	  
	 public static void bhmie()
	 {
	 double amu[]=new  double[100];
	 double theta[]=new  double[100];
	 double pi[]=new  double[100];
	 double tau[]=new  double[100];
	 double pi0[]=new  double[100];
	 double pi1[]=new  double[100];
	 Complex y,xi,xi1,ann,bnn,an,bn,an1,an2,comp1,comp2,comp3,comp4,comp5;
	 Complex d[]=new Complex[3000];
	 double psi0,psi1,psi,dn,dx,xstop,nstop,ymod,dang,rn,chi,chi0,chi1,apsi1,
	 apsi,fn,p,t,qext,qsca,qback,qabs,albedo;
	 int nn,n,nmax,j,jj;
	 dx=xk;
	 an1=new Complex(xk,0.0);
	 ann=new Complex(0.0,0.0);
	 bnn=new Complex(0.0,0.0);
	 an=new Complex(0.0,0.0);
	 bn=new Complex(0.0,0.0);
	 y=an1.times(refrelk);
	 gscatxt=0.0;
	 if(xk<8)
	 {
	 xstop=xk+(4.0)*Math.pow(xk,0.333)+1.0;
	 }
	 else if(xk<=4200)
	 {
	 xstop=xk+(4.5)*Math.pow(xk,0.333)+2.0;
	 }
	 else
	 {
	 xstop=xk+(4.0)*Math.pow(xk,0.333)+2.0;
	 }
	 nstop=xstop;
	 ymod=y.mod();
	 nmax=(int)(Math.max(xstop,ymod))+15;
	 dang=(PI/2.0)/(double)(nangk-1);
	 for(j=1;j<=nangk;j++)
	 {
	 theta[j]=((double)(j)-1.0)*dang;
	 amu[j]=Math.cos(theta[j]);
	 }
	 d[nmax]=new Complex(0.0,0.0);
	 nn=nmax-1;
	 for(n=1;n<=nn;n=n+1)
	 {
	 rn=(double)(nmax-n+1);
	 an1=new Complex(rn,0.0);
	 an2=new Complex(1.0,0.0);
	 d[nmax-n]=an1.div(y).minus(an2.div(d[nmax-n+1].plus(an1.div(y))));
	 }
	 for(j=1;j<=nangk;j++)
	 {
	 pi0[j]=0.0;
	 pi1[j]=1.0;
	 }
	 nn=2*nangk-1;
	 for(j=1;j<=nn;j++)
	 {
	 s1k[j]=new Complex(0.0,0.0);
	 s2k[j]=new Complex(0.0,0.0);
	 }
	 psi0= Math.cos(dx);
	 psi1=Math.sin(dx);
	 chi0= -(Math.sin(xk));
	 chi1=Math.cos(xk);
	 apsi1 =psi1;
	 xi1= new Complex(apsi1,-chi1);
	 qsca=0.0;
	 n=1;
	 boolean flag=true;
	 while(flag)
	 {
	 dn=(double)(n);
	 rn=(double)(n);
	 fn=(2.0*rn+1.0)/(rn*(rn+1.0));
	 psi=(2.0*dn-1.0)*psi1/dx-psi0;
	 apsi=psi;
	 chi=(2.0*rn-1.0)*chi1/(xk)-chi0;
	 xi=new Complex(apsi,-chi);
	 an1=d[n].div(refrelk);
	 comp1=new Complex(rn/(xk),0.0);
	 comp2=new Complex(apsi,0.0);
	 comp3=new Complex(apsi1,0.0);
	 if(n>1)
	 {
	 ann=an;
	 bnn=bn;
	 }
	 an1=an1.plus(comp1);
	 an1=an1.times(comp2);
	 an1=an1.minus(comp3);
	 an2=d[n].div(refrelk);
	 an2=an2.plus(comp1);
	 an2=an2.times(xi);
	 an2=an2.minus(xi1);
	 an=an1.div(an2);
	 an1=refrelk.times(d[n]);
	 an1=an1.plus(comp1);
	 an1=an1.times(comp2);
	 an1=an1.minus(comp3);
	 an2=refrelk.times(d[n]);
	 an2=an2.plus(comp1);
	 an2=an2.times(xi);
	 an2=an2.minus(xi1);
	 bn=an1.div(an2);
	 qsca=qsca+(2.0*rn+1.0)*(Math.pow(an.mod(),2)+Math.pow(bn.mod(),2));
	 gscatxt=gscatxt+(2.0*rn+1.0)*((an.times(bn.conj())).real())/((rn+1.0)*rn);
	 if(n>1)
	 {
	 gscatxt=gscatxt+(rn-1.0)*(rn+1.0)*(ann.times(an.conj()).plus(bnn.times(bn.conj()))).real()/rn;
	 }
	 for(j=1;j<=nangk;j++)
	 {
	 jj=2*nangk-j;
	 pi[j]=pi1[j];
	 tau[j]=rn*amu[j]*pi[j]-(rn+1.0)*pi0[j];
	 p=Math.pow(-1.0,(double)(n-1));
	 comp1=new Complex(pi[j],0.0);
	 comp2=new Complex(tau[j],0.0);
	 comp3=new Complex(fn,0.0);
	 an1=an.times(comp1);
	 an2=bn.times(comp2);
	 an1=comp3.times(an1.plus(an2));
	 s1k[j]=s1k[j].plus(an1);
	 t=Math.pow(-1.0,(double)(n));
	 an1=an.times(comp2);
	 an2=bn.times(comp1);
	 an1=comp3.times(an1.plus(an2));
	 s2k[j]=s2k[j].plus(an1);
	 if(j!=jj)
	 {
	 comp1=new Complex(pi[j],0.0);
	 comp2=new Complex(tau[j],0.0);
	 comp3=new Complex(fn,0.0);
	 comp4=new Complex(p,0.0);
	 comp5=new Complex(t,0.0);
	 s1k[jj]=s1k[jj].plus(comp3.times(an.times(comp1.times(comp4)).plus(bn.times(comp2.times(comp5)))));
	 s2k[jj]=s2k[jj].plus(comp3.times(an.times(comp2.times(comp5)).plus(bn.times(comp1.times(comp4)))));
	 }
	 }
	 psi0=psi1;
	 psi1=psi;
	 apsi1=psi1;
	 chi0=chi1;
	 chi1=chi;
	 xi1=new Complex(apsi1,-chi1);
	 n=n+1;
	 rn=(double)(n);
	 for(j=1;j<=nangk;j++)
	 {
	 pi1[j]=((2.0*rn-1.0)/(rn-1.0))*amu[j]*pi[j];
	 pi1[j]=pi1[j]-rn*pi0[j]/(rn-1.0);
	 pi0[j]=pi[j];
	 }
	 if((n-1-nstop)<0)
	 {
	 flag=true;
	 }
	 else
	 {
	 qsca=(2.0/Math.pow(xk,2))*qsca;
	 qext=(4.0/Math.pow(xk,2))*s1k[1].real();
	 qback=(4.0/Math.pow(xk,2))*Math.pow(s1k[2*(nangk)-1].mod(),2);
	 gscatxt=((4.0/Math.pow(xk, 2))*gscatxt);
	 qscatxt=qsca;
	 qexttxt=qext;
	 qbacktxt=qback;
	 flag=false;
	 }
	 }
	 return;
	 }

	 /*Size distribution functions*/

	 static double dstn(int ntot, double radl, double radh, double rad, double step, int dstn, double rc, 
	 double alfa, double sigma, double rg)
	 {
	 double l,a,b,gama=1.0,nr,rad1,rad2,rad3,rad4,nrn;
	 if(dstn==1)
	 {
	 b=alfa/rc;
	 for(l=alfa;l>=1.0e+0;l=l-1.0e+0)
	 gama= gama*l;
	 a=ntot/(Math.pow(b,-(alfa+1.0))*gama);
	 nr=a*Math.pow(rad,alfa)*Math.exp(-b*rad);
	 }
	 else if(dstn==2)
	 {
	 rad1=rad-rg;
	 rad2=Math.pow(rad1,2.0);
	 rad3=rad2/(2*Math.pow(sigma,2));
	 rad4=1/Math.exp(rad3);
	 nrn=rad4/(Math.pow(2.0*PI,0.5)*sigma);
	 nr=ntot*nrn;
	 }
	 else
	 {
	 rad1=Math.log(rad/rg);
	 rad2=Math.pow(rad1,2.0);
	 rad3=rad2/(2*Math.pow(sigma,2));
	 rad4=1/Math.exp(rad3);
	 nrn=rad4/(rad*(Math.pow(2.0*PI,0.5)*sigma));
	 nr=ntot*nrn;
	 }
	 return(nr);
	 }

	 /*Light scattering calculations on spherical particles*/

	 public static void extra(double stepk,double radlk,double radhk,int ntotk,double rck,double alfak,int dstnk)
	 {
	 double s11tot[]=new double[200];
	 double s12tot[]=new double[200];
	 double s33tot[]=new double[200];
	 double s34tot[]=new double[200];
	 double s11ang[]=new double[200];
	 for(int i=0;i<s11tot.length;i++)
	 {
	 s11tot[i]=0.0;
	 s12tot[i]=0.0;
	 s33tot[i]=0.0;
	 s34tot[i]=0.0;
	 }
	 double sigmak=alfak;
	 double rgk=rck;
	 double ext=0.0;
	 double sca=0.0;
	 double back=0.0;
	 double gr=0.0;
	 double gsc=0.0;
	 xk=0.0;
	 for(double radk=radlk;radk<=radhk;radk=radk+stepk)
	 {
	 double nr=dstn(ntotk,radlk,radhk,radk,stepk,dstnk,rck,alfak,sigmak,rgk);
	 nr=stepk*nr;
	 xk=2.0*PI*radk*(refmed.real())/wavel;
	 dang=(PI/2.0)/(double)(nangk-1);
	 bhmie();
	 gr=gr+PI*Math.pow(radk,2)*nr;
	 sca=sca+PI*Math.pow(radk,2)*nr*qscatxt;
	 ext=ext+PI*Math.pow(radk,2)*nr*qexttxt;
	 back=back+PI*Math.pow(radk,2)*nr*qbacktxt;
	 gsc=gsc+PI*Math.pow(radk,2)*nr*gscatxt;
	 qscatxt=sca/gr;
	 qexttxt=ext/gr;
	 qbacktxt=back/gr;
	 gscatxt=gsc/gr;
	 gtxt=gscatxt/qscatxt;
	 qprtxt=qexttxt-gscatxt;
	 qabstxt=qexttxt-qscatxt;
	 albedotxt=qscatxt/qexttxt;
	 nan=2*nangk-1;
	 for(j=1;j<=nan;j++)
	 {
	 s11=0.5e+0*(Math.pow(s2k[j].mod(),2.0)+Math.pow(s1k[j].mod(),2.0));
	 s11tot[j]=s11tot[j]+s11*nr;
	 s12=(0.5e+0)*(Math.pow(s2k[j].mod(),2.0)-Math.pow((s1k[j].mod()),2.0));
	 s12tot[j]=s12tot[j]+s12*nr;
	 s33=(s2k[j].times(s1k[j].conj())).real();
	 s33tot[j]=s33tot[j]+s33*nr;
	 s34=(s2k[j].times((s1k[j].conj()))).imag();
	 s34tot[j]=s34tot[j]+s34*nr;
	 }
	 }
	 int count=0;
	 for(j=1;j<=nan;j++)
	 {
	 ang=dang*((double)j-1.0e+0)*(180.0e+0/PI);
	 ss11[count]=s11tot[j]/s11tot[1];
	 ss12[count]=-s12tot[j]/s11tot[j];
	 ss33[count]=s33tot[j]/s11tot[j];
	 ss34[count]=s34tot[j]/s11tot[j];
	 count++;
	 }
	
	 }

	

	public String getMessage() {
		return message;
	}

	public void setMessage(String message) {
		this.message = message;
	}




	public CartesianChartModel getLinearModel2() {
		return linearModel2;
	}




	public CartesianChartModel getLinearModel3() {
		return linearModel3;
	}




	public CartesianChartModel getLinearModel4() {
		return linearModel4;
	}
	
	public CartesianChartModel getLinearModel1() {
		return linearModel1;
	}




	public String getFixedrrefractiveIndexTextBox() {
		return fixedrrefractiveIndexTextBox;
	}




	public void setFixedrrefractiveIndexTextBox(String fixedrrefractiveIndexTextBox) {
		this.fixedrrefractiveIndexTextBox = fixedrrefractiveIndexTextBox;
	}




	public String getFixedirefractiveIndexTextBox() {
		return fixedirefractiveIndexTextBox;
	}




	public void setFixedirefractiveIndexTextBox(String fixedirefractiveIndexTextBox) {
		this.fixedirefractiveIndexTextBox = fixedirefractiveIndexTextBox;
	}




	public boolean isRealindexRadioButton() {
		return realindexRadioButton;
	}




	public void setRealindexRadioButton(boolean realindexRadioButton) {
		this.realindexRadioButton = realindexRadioButton;
	}




	public boolean isImagindexRadioButton() {
		return imagindexRadioButton;
	}




	public void setImagindexRadioButton(boolean imagindexRadioButton) {
		this.imagindexRadioButton = imagindexRadioButton;
	}




	public boolean isSphereRadioButton() {
		return sphereRadioButton;
	}




	public void setSphereRadioButton(boolean sphereRadioButton) {
		this.sphereRadioButton = sphereRadioButton;
	}




	public boolean isSpheroidRadioButton() {
		return spheroidRadioButton;
	}




	public void setSpheroidRadioButton(boolean spheroidRadioButton) {
		this.spheroidRadioButton = spheroidRadioButton;
	}




	public boolean isCyclinderRadioButton() {
		return cyclinderRadioButton;
	}




	public void setCyclinderRadioButton(boolean cyclinderRadioButton) {
		this.cyclinderRadioButton = cyclinderRadioButton;
	}




	public String getSelectsizeRef() {
		return selectsizeRef;
	}




	public void setSelectsizeRef(String selectsizeRef) {
		this.selectsizeRef = selectsizeRef;
	}





	public boolean isBoolsizeRI() {
		return boolsizeRI;
	}





	public void setBoolsizeRI(boolean boolsizeRI) {
		this.boolsizeRI = boolsizeRI;
	}





	public String getMinRadiusTextBox() {
		return minRadiusTextBox;
	}





	public void setMinRadiusTextBox(String minRadiusTextBox) {
		this.minRadiusTextBox = minRadiusTextBox;
	}





	public String getMaxRadiusTextBox() {
		return maxRadiusTextBox;
	}





	public void setMaxRadiusTextBox(String maxRadiusTextBox) {
		this.maxRadiusTextBox = maxRadiusTextBox;
	}





	public String getStepRadiusTextBox() {
		return stepRadiusTextBox;
	}





	public void setStepRadiusTextBox(String stepRadiusTextBox) {
		this.stepRadiusTextBox = stepRadiusTextBox;
	}





	public boolean isBoolsizeRISpCy() {
		return boolsizeRISpCy;
	}





	public void setBoolsizeRISpCy(boolean boolsizeRISpCy) {
		this.boolsizeRISpCy = boolsizeRISpCy;
	}

	public String getSelectcomboTypesPara() {
		return selectcomboTypesPara;
	}

	public void setSelectcomboTypesPara(String selectcomboTypesPara) {
		this.selectcomboTypesPara = selectcomboTypesPara;
	}

	public boolean isBoolsizeABCD() {
		return boolsizeABCD;
	}

	public void setBoolsizeABCD(boolean boolsizeABCD) {
		this.boolsizeABCD = boolsizeABCD;
	}

	public String getMinAccuracyTextBox() {
		return minAccuracyTextBox;
	}

	public void setMinAccuracyTextBox(String minAccuracyTextBox) {
		this.minAccuracyTextBox = minAccuracyTextBox;
	}

	public String getMaxAccuracyTextBox() {
		return maxAccuracyTextBox;
	}

	public void setMaxAccuracyTextBox(String maxAccuracyTextBox) {
		this.maxAccuracyTextBox = maxAccuracyTextBox;
	}

	public String getStepAccuracyTextBox() {
		return stepAccuracyTextBox;
	}

	public void setStepAccuracyTextBox(String stepAccuracyTextBox) {
		this.stepAccuracyTextBox = stepAccuracyTextBox;
	}

	public String getEqvRadiusTextBox() {
		return eqvRadiusTextBox;
	}

	public void setEqvRadiusTextBox(String eqvRadiusTextBox) {
		this.eqvRadiusTextBox = eqvRadiusTextBox;
	}

	public String getDisplayTextABVE() {
		return displayTextABVE;
	}

	public void setDisplayTextABVE(String displayTextABVE) {
		this.displayTextABVE = displayTextABVE;
	}

	public String getDisplayTextComBo() {
		return displayTextComBo;
	}

	public void setDisplayTextComBo(String displayTextComBo) {
		this.displayTextComBo = displayTextComBo;
	}

	public boolean isBoolrefractive() {
		return boolrefractive;
	}

	public void setBoolrefractive(boolean boolrefractive) {
		this.boolrefractive = boolrefractive;
	}

	public String getDisplayTextRealImag() {
		return displayTextRealImag;
	}

	public void setDisplayTextRealImag(String displayTextRealImag) {
		this.displayTextRealImag = displayTextRealImag;
	}

	public String getDisplayTextLegen() {
		return displayTextLegen;
	}

	public void setDisplayTextLegen(String displayTextLegen) {
		this.displayTextLegen = displayTextLegen;
	}
	public String getSelectcomboTypesShapImg() {
		return selectcomboTypesShapImg;
	}
	public void setSelectcomboTypesShapImg(String selectcomboTypesShapImg) {
		this.selectcomboTypesShapImg = selectcomboTypesShapImg;
	}
	public String getRadiusImgRealTextBox() {
		return radiusImgRealTextBox;
	}
	public void setRadiusImgRealTextBox(String radiusImgRealTextBox) {
		this.radiusImgRealTextBox = radiusImgRealTextBox;
	}
	public String getMinirefractiveIndexTextBox() {
		return minirefractiveIndexTextBox;
	}
	public void setMinirefractiveIndexTextBox(String minirefractiveIndexTextBox) {
		this.minirefractiveIndexTextBox = minirefractiveIndexTextBox;
	}
	public String getMaxrefractiveIndexTextBox() {
		return maxrefractiveIndexTextBox;
	}
	public void setMaxrefractiveIndexTextBox(String maxrefractiveIndexTextBox) {
		this.maxrefractiveIndexTextBox = maxrefractiveIndexTextBox;
	}
	public String getSteprefractiveIndexTextBox() {
		return steprefractiveIndexTextBox;
	}
	public void setSteprefractiveIndexTextBox(String steprefractiveIndexTextBox) {
		this.steprefractiveIndexTextBox = steprefractiveIndexTextBox;
	}
	public boolean isBoolrefPara() {
		return boolrefPara;
	}
	public void setBoolrefPara(boolean boolrefPara) {
		this.boolrefPara = boolrefPara;
	}
	public String getDisplayTextLast() {
		return displayTextLast;
	}
	public void setDisplayTextLast(String displayTextLast) {
		this.displayTextLast = displayTextLast;
	}
	public String getDisplayexpectedRange() {
		return displayexpectedRange;
	}
	public void setDisplayexpectedRange(String displayexpectedRange) {
		this.displayexpectedRange = displayexpectedRange;
	}



	public StreamedContent getFile() {
		return file;
	}



	public void setFile(StreamedContent file) {
		this.file = file;
	}



	public List<TuScatPojo> getListofdata() {
		return listofdata;
	}



	public void setListofdata(List<TuScatPojo> listofdata) {
		this.listofdata = listofdata;
	}



	public String getMinimunError() {
		return minimunError;
	}



	public void setMinimunError(String minimunError) {
		this.minimunError = minimunError;
	}
	public boolean isBoolsizeRISproid() {
		return boolsizeRISproid;
	}
	public void setBoolsizeRISproid(boolean boolsizeRISproid) {
		this.boolsizeRISproid = boolsizeRISproid;
	}
	public String getDisplayTextComBoEqui() {
		return displayTextComBoEqui;
	}
	public void setDisplayTextComBoEqui(String displayTextComBoEqui) {
		this.displayTextComBoEqui = displayTextComBoEqui;
	}
	public boolean isBoolsizeABCD1() {
		return boolsizeABCD1;
	}
	public void setBoolsizeABCD1(boolean boolsizeABCD1) {
		this.boolsizeABCD1 = boolsizeABCD1;
	}
	public boolean isBoolrefrctiveSphere() {
		return boolrefrctiveSphere;
	}
	public void setBoolrefrctiveSphere(boolean boolrefrctiveSphere) {
		this.boolrefrctiveSphere = boolrefrctiveSphere;
	}
	public boolean isBoolrefrctiveCylinder() {
		return boolrefrctiveCylinder;
	}
	public void setBoolrefrctiveCylinder(boolean boolrefrctiveCylinder) {
		this.boolrefrctiveCylinder = boolrefrctiveCylinder;
	}
	public boolean isBoolrefrctiveSpheroid() {
		return boolrefrctiveSpheroid;
	}
	public void setBoolrefrctiveSpheroid(boolean boolrefrctiveSpheroid) {
		this.boolrefrctiveSpheroid = boolrefrctiveSpheroid;
	}
	
	public String getSelectParameter() {
		return selectParameter;
	}
	public void setSelectParameter(String selectParameter) {
		this.selectParameter = selectParameter;
	}
	public boolean isBoolsizeABCDspheroidab() {
		return boolsizeABCDspheroidab;
	}
	public void setBoolsizeABCDspheroidab(boolean boolsizeABCDspheroidab) {
		this.boolsizeABCDspheroidab = boolsizeABCDspheroidab;
	}
	public boolean isBoolsizeABCDspheroiddl() {
		return boolsizeABCDspheroiddl;
	}
	public void setBoolsizeABCDspheroiddl(boolean boolsizeABCDspheroiddl) {
		this.boolsizeABCDspheroiddl = boolsizeABCDspheroiddl;
	}
	public String getSelectparameter2() {
		return selectparameter2;
	}
	public void setSelectparameter2(String selectparameter2) {
		this.selectparameter2 = selectparameter2;
	}



	public FacesMessage getMassage() {
		return massage;
	}



	public void setMassage(FacesMessage massage) {
		this.massage = massage;
	}





}
