package coreservlets;

import java.util.HashMap;

public class MinimumErrorPojo {
	
	private HashMap<String, String> s11;
	private HashMap<String, String> s12;
	private HashMap<String, String> s33;
	private HashMap<String, String> s34;
	
	public MinimumErrorPojo()
	{
		s11=new HashMap<String, String>();
		s12=new HashMap<String, String>();
		s33=new HashMap<String, String>();
		s34=new HashMap<String, String>();
	}
	
	
	public  void addInformation(String seq,String ss11,String ss12,String ss33,String ss34)
	{
		if(null!=seq&&null!=ss11&&null!=ss11&&null!=ss33&&null!=ss34)
		{
			s11.put(seq.trim(),ss11.trim());
			s12.put(seq.trim(),ss12.trim());
			s33.put(seq.trim(),ss33.trim());
			s34.put(seq.trim(),ss34.trim());
		}
	}
	
	
	
	
	
	
	
	
	
	public HashMap<String, String> getS11() {
		return s11;
	}
	public void setS11(HashMap<String, String> s11) {
		this.s11 = s11;
	}
	public HashMap<String, String> getS12() {
		return s12;
	}
	public void setS12(HashMap<String, String> s12) {
		this.s12 = s12;
	}
	public HashMap<String, String> getS33() {
		return s33;
	}
	public void setS33(HashMap<String, String> s33) {
		this.s33 = s33;
	}
	public HashMap<String, String> getS34() {
		return s34;
	}
	public void setS34(HashMap<String, String> s34) {
		this.s34 = s34;
	}

}
