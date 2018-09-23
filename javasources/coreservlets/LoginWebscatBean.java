package coreservlets;

import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import javax.faces.application.FacesMessage;
import javax.faces.bean.ManagedBean;
import javax.faces.bean.SessionScoped;
import javax.faces.context.FacesContext;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;

import org.primefaces.context.RequestContext;

@ManagedBean(name="LoginWebscatBean")
@SessionScoped
public class LoginWebscatBean {
	
	
	private String emailAddress;
 	private String passwordWebscat;
 	
 /*	public LoginWebscatBean()
 	{
 		
 	}*/
 	
 	
 	
 	/// Main login script here
 	public void loginwebsactapplication() 
 	{
 		try
 		{
 			Class.forName("com.mysql.jdbc.Driver").newInstance();
			//String url="jdbc:mysql://127.11.145.130:3306/webscat";
 			String url = "jdbc:mysql://localhost:3306/webscat";
 			Connection conn = DriverManager.getConnection(url,
 					"webscat", "webscat");
 		PreparedStatement pstmt = null;
 		ResultSet rs = null;
 		
 			if(conn != null)
 			{
 				
 				pstmt = conn.prepareStatement("select EMAILADDRESS ,CPASSWORD from LOGINTABLE");			
 				rs = pstmt.executeQuery();
 				while (rs.next()) 
 				{
 					System.out.println("Hello ! I m working"+rs.getString("EMAILADDRESS")+"\n"+rs.getString("cpassword"));
 					
 					String emaildb=rs.getString("EMAILADDRESS");
 					String passworddb=rs.getString("cpassword");
 					if(null != emailAddress&&emaildb.equalsIgnoreCase(emailAddress))
 					{
 						if(null != passwordWebscat&&passworddb.equalsIgnoreCase(passwordWebscat))
 						{
 							//FacesContext context = FacesContext.getCurrentInstance();
 							//context = FacesContext.getCurrentInstance();
 							//context.getExternalContext().redirect(context.getExternalContext().getRequestContextPath() + "/aboutus.jsf");
 					       
 							FacesContext context1 = FacesContext.getCurrentInstance();
 							context1.getExternalContext().redirect(context1.getExternalContext().getRequestContextPath() + "/aboutus.jsf");
 						    HttpServletRequest reques = (HttpServletRequest)context1.getExternalContext().getRequest();  
 						    HttpSession session = reques.getSession(true);
 							session.setAttribute("emailid",emailAddress);
 							String email = (String) session.getAttribute("emailid");
 							
 							
 							/*RequestContext context = RequestContext.getCurrentInstance();  
 	 						FacesMessage msg22 = null;
 	 						msg22 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "Successfull login"+"..."+email);
 	 						FacesContext.getCurrentInstance().addMessage(null, msg22);*/
 						}
 						else if(null == passwordWebscat)
 						{
 							RequestContext context = RequestContext.getCurrentInstance();  
 	 						FacesMessage msg3 = null;
 	 						msg3 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "Password can't be empty");
 	 						FacesContext.getCurrentInstance().addMessage(null, msg3);
 						}
 						else if(!passworddb.equalsIgnoreCase(passwordWebscat))
 						{
 							RequestContext context = RequestContext.getCurrentInstance();  
 	 						FacesMessage msg5 = null;
 	 						msg5 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "Wrong password");
 	 						FacesContext.getCurrentInstance().addMessage(null, msg5);
 						}
 					}
 					else if(null == emailAddress)
 					{
 						RequestContext context = RequestContext.getCurrentInstance();  
 						FacesMessage msg2 = null;
 						msg2 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "User name can't be empty");
 						FacesContext.getCurrentInstance().addMessage(null, msg2);
 					}
 					else if(!emaildb.equalsIgnoreCase(emailAddress))
 					{
 						RequestContext context = RequestContext.getCurrentInstance();  
 						FacesMessage msg4 = null;
 						msg4 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "Wrong user name");
 						FacesContext.getCurrentInstance().addMessage(null, msg4);
 					}
 				}
 		
 		
 			}
 			else if(conn==null)
 			{
 				System.out.println("Can not connect to the database .....");
 				RequestContext context = RequestContext.getCurrentInstance();  
				FacesMessage msg1 = null;
				msg1 = new FacesMessage(FacesMessage.SEVERITY_INFO, "Server Says:", "Can't connect to the server");
				FacesContext.getCurrentInstance().addMessage(null, msg1);
 			}
 		}
 		catch(Exception e)
 		{
 			e.printStackTrace();
 			System.out.println("e.printStackTrace()");
 		}
 		
 		
 	}
 	
 	
	public String getEmailAddress() {
		return emailAddress;
	}
	public void setEmailAddress(String emailAddress) {
		this.emailAddress = emailAddress;
	}
	public String getPasswordWebscat() {
		return passwordWebscat;
	}
	public void setPasswordWebscat(String passwordWebscat) {
		this.passwordWebscat = passwordWebscat;
	}

}
