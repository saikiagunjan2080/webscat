package util;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;

import javax.faces.context.FacesContext;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpSession;
public class Diff {
	public static void main(String[] args) throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException, IOException {


		
		Connection con=null;
		Class.forName("com.mysql.jdbc.Driver").newInstance();
		//String url="jdbc:mysql:/$OPENSHIFT_MYSQL_DB_HOST:$OPENSHIFT_MYSQL_DB_PORT/webscat";
		String url="jdbc:mysql://localhost:3306/gyanbharatidb";
		//String url = "jdbc:mysql://localhost:3306/sciencedb";
		Connection conn = DriverManager.getConnection(url,"gyanbharatidb", "gyanbharatidb");
		PreparedStatement pstmt = null;
 		ResultSet rs = null;
 		/*if(conn!=null)
 		{
 			pstmt = conn.prepareStatement("select EMAILADDRESS ,CPASSWORD from LOGINTABLE");			
			rs = pstmt.executeQuery();
			while (rs.next())
			{
			String emaildb=rs.getString("EMAILADDRESS");
			String passworddb=rs.getString("cpassword");
			System.out.println("add"+emaildb+"....."+passworddb);
			}
			try
			{
				//FacesContext context1 = FacesContext.getCurrentInstance();
				//context1.getExternalContext().redirect(context1.getExternalContext().getRequestContextPath() + "http://webscat-tewebapplication.rhcloud.com/WEBSCATwebapplication/aboutus.jsf");	
			}
			catch(Exception e)
	 		{
	 			e.printStackTrace();
	 		}
			
			    
 		}*/
		

		System.out.println("connected");
		}

}
