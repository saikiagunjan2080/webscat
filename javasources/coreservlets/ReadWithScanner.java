package coreservlets;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.Scanner;



public class ReadWithScanner {
	
	MinimumErrorPojo minimumErrorPojo;
	
	public void addAllPoint()
	{
	      String ss11 = "";
	      String ss12 = "";
	      String ss33 = "";
	      String ss34 = "";
	      
		for (int i=0;i<181;i++) {
			String seq=""+i;
			ss11=this.minimumErrorPojo.getS11().get(seq.trim());
			ss12=this.minimumErrorPojo.getS12().get(seq.trim());
			ss33=this.minimumErrorPojo.getS33().get(seq.trim());
			ss34=this.minimumErrorPojo.getS34().get(seq.trim());
			//System.out.println(seq+" |,| "+ ss11+" |,| "+ ss12+" |,| "+ss33+" |,| "+ss34);
		}
	}
	
	public void addAllPoint(String file)
	{
		try{
			FileInputStream fstream = new FileInputStream(file);
			BufferedReader br = new BufferedReader(new InputStreamReader(fstream));
			ReadWithScanner readWithScanner=new ReadWithScanner();;
			String strLine;
			//Read File Line By Line
			while ((strLine = br.readLine()) != null)   {
			  // Print the content on the console
			  //System.out.println (strLine);
				readWithScanner.processLine(strLine);
			}
			readWithScanner.addAllPoint();
			//Close the input stream
			br.close();
			}
			catch (Exception e) {
				// TODO: handle exception
				e.printStackTrace();
			}
	}
	
	
	public static void main(String args[])
	{
		try{
		FileInputStream fstream = new FileInputStream("C:\\Users\\pritom\\Downloads\\Datatobeinserted.txt");
		BufferedReader br = new BufferedReader(new InputStreamReader(fstream));
		ReadWithScanner readWithScanner=new ReadWithScanner();;
		String strLine;
		//Read File Line By Line
		while ((strLine = br.readLine()) != null)   {
		  // Print the content on the console
		  //System.out.println (strLine);
			readWithScanner.processLine(strLine);
		}
		readWithScanner.addAllPoint();
		//Close the input stream
		br.close();
		}
		catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
		}
	}
	
	
	
	
	public ReadWithScanner()
	{
		this.minimumErrorPojo=new MinimumErrorPojo();
	}
	public final MinimumErrorPojo processLineByLine(String topology)
	{
		//Note that FileReader is used, not File, since File is not Closeable
	    Scanner scanner = new Scanner(topology);
	    try {
	      //first use a Scanner to get each line
	      while ( scanner.hasNextLine() ){
	    	 //System.out.println("test:processLineByLine");
	        processLineByLine( scanner.nextLine());
	      }
	    }
	    finally {
	      //ensure the underlying stream is always closed
	      //this only has any effect if the item passed to the Scanner
	      //constructor implements Closeable (which it does in this case).
	      scanner.close();
	    }
	    return minimumErrorPojo;
	}
	
	protected void processLine(String aLine){
	    //use a second Scanner to parse the content of each line 
	    Scanner scanner = new Scanner(aLine);
	   
	    scanner.useDelimiter(",");
	    if ( scanner.hasNext() ){
	      String seq = scanner.next();
	      String ss11 = scanner.next();
	      String ss12 = scanner.next();
	      String ss33 = scanner.next();
	      String ss34 = scanner.next();
	     
	      this.minimumErrorPojo.addInformation(seq, ss11, ss12,ss33,ss34);
	      //System.out.println(seq+" |,| "+ ss11+" |,| "+ ss12+" |,| "+ss33+" |,| "+ss34);
				
	    }
	    else {
	      log("Empty or invalid line. Unable to process.");
	    }
	  
	    //no need to call scanner.close(), since the source is a String
	  }
	private void log(Object aObject) {
		//System.out.println(String.valueOf(aObject));
		
	}
}
