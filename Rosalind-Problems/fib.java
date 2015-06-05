import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.math.BigInteger;

/*
 *    (..)If F_n represents the number of rabbit pairs alive after the n-th month,
 *    then we obtain the Fibonacci sequence having terms F_n that are defined
 *    by the recurrence relation F_n=F_n−1+F_n−2 (with F1=F2=1 to initiate the
 *    sequence). Although the sequence bears Fibonacci's name, it was known 
 *    to Indian mathematicians over two millennia ago.(..)
	
	Given: Positive integers n≤40 and k≤5.
	Return: The total number of rabbit pairs that will be present after n months
	if each pair of reproduction-age rabbits produces a litter of k rabbit pairs
	in each generation (instead of only 1 pair).
*/


public class fib {
  
	
	public static BigInteger getBInt(String val){
		return new BigInteger(val);
		}
	
	public static BigInteger[] Reader(){
		File file = new File("datasets/rosalind_fib.txt");
		String pattern = "(\\d+) (\\d)";
				
		try {
			Scanner in = new Scanner(file);
			String line = in.nextLine();

			BigInteger n = getBInt(line.replaceAll(pattern,"$1"));
			BigInteger k = getBInt(line.replaceAll(pattern,"$2"));
		    
		  BigInteger[] a = {n,k}; 
			return a;
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return null; 
	    
		
	}
	
  
  
    public static BigInteger getFib(BigInteger n, BigInteger k) {
    	if(n.intValue()==0) return getBInt("0");
    	else if (n.intValue() <= 2) return getBInt("1");
        else return getFib(n.subtract(getBInt("1")), k).add((k.multiply(getFib(n.subtract(getBInt("2")), k))));
    }



    
    public static void main(String[] args) {
    	
    	BigInteger[] r = Reader();
    	BigInteger n = r[0];
    	BigInteger k = r[1];
    	
      System.out.println(getFib(n,k));
      
    }
}

