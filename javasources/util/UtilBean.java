package util;


import java.util.Scanner;


public class UtilBean {

	public static void main(String[] args)
	{
		int num1,num2,sum;
		Scanner in=new Scanner(System.in);
		System.out.println("Please give Value of number1\n");
		num1=in.nextInt();
		System.out.println("Please give Value of number2\n");
		num2=in.nextInt();
		sum=num1+num2;
		System.out.println("Result::"+sum);
		
	}
	
}
