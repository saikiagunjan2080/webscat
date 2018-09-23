package coreservlets;

public class TuScatPojo {
	
	private double angle;
	private double s11P;
	private double s12P;
	private double s33P;
	private double s34P;
	
	
	public TuScatPojo(double angle,double s11P,double s12P,double s33P,double s34P)
	{
		this.angle=angle;
		this.s11P=s11P;
		this.s12P=s12P;
		this.s33P=s33P;
		this.s34P=s34P;
		
	}


	public double getAngle() {
		return angle;
	}


	public void setAngle(double angle) {
		this.angle = angle;
	}


	public double getS11P() {
		return s11P;
	}


	public void setS11P(double s11p) {
		s11P = s11p;
	}


	public double getS12P() {
		return s12P;
	}


	public void setS12P(double s12p) {
		s12P = s12p;
	}


	public double getS33P() {
		return s33P;
	}


	public void setS33P(double s33p) {
		s33P = s33p;
	}


	public double getS34P() {
		return s34P;
	}


	public void setS34P(double s34p) {
		s34P = s34p;
	}


	
	

}
