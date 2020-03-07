import java.awt.List;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;
import java.util.Vector;
import java.util.stream.Collectors;

//import ConcensusTest.AminoAcid;

public class MotifExtraction1 {
	public static enum AminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, Dash};
	public static String aminoAcidLets = "ARNDCQEGHILKMFPSTWYV";
	
	public static void main(String [] args ) throws FileNotFoundException {
		long start = System.currentTimeMillis();
		File file = new File("Motif1.txt");
		Scanner sc = new Scanner (file);
		String input = sc.next();
		Vector<String> Peptides = new Vector<String>();
		for(int i = 0; i<input.length(); i+=7) {
			Peptides.add(input.substring(i,i+7));
		}
		for(int m=0; m < Peptides.size(); m++) System.out.println(m+":"+ " " + Peptides.get(m));
		ArrayList<Set<Integer>> AdjacencyList = new ArrayList();
		AdjacencyList = AdjList(Peptides);
		
		ArrayList <Set<Integer>> Result = new ArrayList();
    	Result = (ArrayList<Set<Integer>>) findCliques(AdjacencyList).stream().filter(c -> c.size() > 1).collect(Collectors.toList());
    	
    	System.out.println();
    	System.out.println("We have:  " + Result.size() + " "+ " Cliques in the Graph");
    	System.out.println();
    	
    	// Writing to a file
    	
    	File file2 = new File("Result1.txt");
		PrintStream out = new PrintStream(file2);
		out.println();
		out.println();
		out.println("          "+ "Predefined Motif is: " + "A--TG-Y");
	//	out.println("          "+ "Predefined Motif is: " + "MR--E-Y");
	//	out.println("          "+ "Predefined Motif is: " + "Q-W--HN");
	//	out.println("          "+ "Predefined Motif is: " + "---LSTM");
	//	out.println("          "+ "Predefined Motif is: " + "-R-CI-V");
	//	out.println("          "+ "Predefined Motif is: " + "MK-F--P");
	//	out.println("          "+ "Predefined Motif is: " + "A-D-KT-");


		out.println();
		out.println("          "+ "Generating Peptides (4Mers) :");
		out.println();
		for(int m=0; m < Peptides.size(); m++) out.println(m+":"+ " " + Peptides.get(m));
		out.println();
		
    	out.println("We have:  " + Result.size() + " "+ " Cliques in the Graph");
    	out.println();
    	
    	// Target Motifs
    	ArrayList<AminoAcid> TM1 = new ArrayList<AminoAcid>();
    //	ArrayList<AminoAcid> TM2 = new ArrayList<AminoAcid>();
    //	ArrayList<AminoAcid> TM3 = new ArrayList<AminoAcid>();
    //	ArrayList<AminoAcid> TM4 = new ArrayList<AminoAcid>();
    //	ArrayList<AminoAcid> TM5 = new ArrayList<AminoAcid>();
    //	ArrayList<AminoAcid> TM6 = new ArrayList<AminoAcid>();
    //	ArrayList<AminoAcid> TM7 = new ArrayList<AminoAcid>();

    	TM1.add(AminoAcid.A);TM1.add(AminoAcid.Dash);TM1.add(AminoAcid.Dash);TM1.add(AminoAcid.T);TM1.add(AminoAcid.G);TM1.add(AminoAcid.Dash);TM1.add(AminoAcid.Y);
    //	TM2.add(AminoAcid.M);TM2.add(AminoAcid.R);TM2.add(AminoAcid.Dash);TM2.add(AminoAcid.Dash);TM2.add(AminoAcid.E);TM2.add(AminoAcid.Dash);TM2.add(AminoAcid.Y);
    //	TM3.add(AminoAcid.Q);TM3.add(AminoAcid.Dash);TM3.add(AminoAcid.W);TM3.add(AminoAcid.Dash);TM3.add(AminoAcid.Dash);TM3.add(AminoAcid.H);TM3.add(AminoAcid.N);
    //	TM4.add(AminoAcid.Dash);TM4.add(AminoAcid.Dash);TM4.add(AminoAcid.Dash);TM4.add(AminoAcid.L);TM4.add(AminoAcid.S);TM4.add(AminoAcid.T);TM4.add(AminoAcid.M);
    //	TM5.add(AminoAcid.Dash);TM5.add(AminoAcid.R);TM5.add(AminoAcid.Dash);TM5.add(AminoAcid.C);TM5.add(AminoAcid.I);TM5.add(AminoAcid.Dash);TM5.add(AminoAcid.V);
    //	TM6.add(AminoAcid.M);TM6.add(AminoAcid.K);TM6.add(AminoAcid.Dash);TM6.add(AminoAcid.F);TM6.add(AminoAcid.Dash);TM6.add(AminoAcid.Dash);TM6.add(AminoAcid.P);
    //	TM6.add(AminoAcid.A);TM7.add(AminoAcid.Dash);TM7.add(AminoAcid.D);TM7.add(AminoAcid.Dash);TM7.add(AminoAcid.K);TM7.add(AminoAcid.T);TM7.add(AminoAcid.Dash);
    	
    	int Cnt = 1;
    	int nMotifs = 1;
    	int CntNew = CountNum(Cnt, Result);
    	
    	int M1 = 0; int M2 = 0; int M3 = 0; int M4 = 0; int M5 = 0;int M6 = 0;int M7 = 0;
    	int C = 0;
    	int falsePositive = 0;
    	for(int j = 0; j < CntNew; j++) {
    	
    		System.out.println("      Round: " + " " + j);
    		out.println("      Round: " + " " + j);
    		System.out.println();
    		out.println();
    		
    		int MaxCliqueSize = 0;
    		// k = Index of Maximal Clique in the ArrayList of All Cliques(Result)
    		int k = 0;
    		for(int i = 0; i < Result.size(); i++) {
    			if (Result.get(i).size() > MaxCliqueSize) { 
    				MaxCliqueSize = Result.get(i).size();
    				k = i;
    			}
    		}   	
    		System.out.println("Size of the Maximum Clique in the Graph is " + " "  + Result.get(k).size());	
    		out.println("Size of the Maximum Clique in the Graph is " + " "  + Result.get(k).size());

    		ArrayList<String> Final4Mers = new ArrayList<String>();
    			for(Iterator it = Result.get(k).iterator(); it.hasNext();) {
    	    		Final4Mers.add(Peptides.get((Integer)it.next()));    
    	    	}
    				System.out.println("Final Answer is" + ":" );
    				out.println("Final Answer is" + ":" );
    				System.out.println();
    				out.println();
    				System.out.println("    " + Consensus(Final4Mers).toString());
    				out.println("    " + Consensus(Final4Mers).toString());
    				Consensus(Final4Mers);
    				
    				if(Consensus(Final4Mers).get(0) == AminoAcid.A && Consensus(Final4Mers).get(3) == AminoAcid.T && 
    						Consensus(Final4Mers).get(4) == AminoAcid.G && Consensus(Final4Mers).get(6) == AminoAcid.Y) {
    					M1+=1;
    				}
    				else {
    					falsePositive += 1;
    				}
    			/*	if(Consensus(Final4Mers).equals(TM3)) {
    					M3+=1;
    				}
    			*/	
    				
    				
    				Result.remove(k);
    				System.out.println();
    				System.out.println();
    				out.println();
    	}
    	//while
    	
		
		if(M1 > 1) M1 = 1;
		if(M2 > 1) M2 = 1;
		if(M3 > 1) M3 = 1;
		if(M4 > 1) M4 = 1;
		if(M5 > 1) M5 = 1;
		if(M6 > 1) M6 = 1;
		if(M7 > 1) M7 = 1;
		
		C = M1;
		
    	double Specificity = C*100/ (nMotifs+falsePositive);
		double Sensitivity = C*100/ nMotifs;
		
		System.out.println();
		System.out.println("     Number of Cliques is: " + " " + CntNew);
		System.out.println("     Sensitivity is: " + "%" + Sensitivity);
		System.out.println("     Specificity is: " + "%" + Specificity);
		System.out.println("     -------------------------");
		
		out.println();
		out.println("     Number of Cliques is: " + " " + CntNew);
		out.println("     Sensitivity is: " + "%" + Sensitivity);
		out.println("     Specificity is: " + "%" + Specificity);
		out.println("     -------------------------");
		long finish = System.currentTimeMillis();
		long timeElapsed = finish - start;
		System.out.println("Time Elapsed: " + " " + timeElapsed);
	}
	
	public static ArrayList<Set<Integer>> AdjList(Vector<String> seq){
		ArrayList<Set<Integer>> adjacencyList = new ArrayList<Set<Integer>> ();
		for(int i = 0; i < seq.size(); i++) {
			adjacencyList.add(new HashSet());
		}
		for(int i = 0; i < seq.size(); i++) 
		{
			for(int j = 0; j < seq.size(); j++) 
			{
				if(i == j)continue;
				int score = 0;
				for(int k = 0; k < 7; k++) 
				{
					if( seq.get(i).charAt(k) == seq.get(j).charAt(k) && seq.get(i).charAt(k) != '-') 
					{
						score+=1;
					}
					
				}
				if (score >= 3) 
				{
					adjacencyList.get(i).add(j);
				//	adjacencyList.get(j).add(i);
				}
				
			//	if(adjacencyList.get(i).contains(j))adjacencyList.get(j).add(i);
			}
		}
		
		return adjacencyList;

	}

	public static int[] Repetition(Vector<String> vect) {
		int[] frequency = new int[vect.size()];
		for (int k = 0; k < frequency.length; k++) {
			frequency[k] = 1;
		}
		for(int i = 0; i < frequency.length-1; i++) {
			for (int j = i+1; j < frequency.length; j++) {
				if (vect.get(i).equals(vect.get(j))){
					frequency[i] += 1;
				}
			}
		}
		return frequency;
	}

	/**
	 * Simple implementation of Bron - Kerbosch algorithm for finding all maximum cliques (psedocode is on Wiki)
	 */
	public static Set<Set<Integer>> findCliques(ArrayList<Set<Integer>> adjacencyList) {
	    Set<Integer> p = new HashSet<>();
	    for (int i = 0; i < adjacencyList.size(); i++) {
	        p.add(i);
	    }
	    return findCliques(new HashSet<>(), p, new HashSet<>(), adjacencyList);
	}

	private static Set<Set<Integer>> findCliques(Set<Integer> clique, Set<Integer> p, Set<Integer> x, ArrayList<Set<Integer>> adjacencyList) {
	    Set<Set<Integer>> result = new HashSet<>();
	    if (p.isEmpty() && x.isEmpty()) {
	        result.add(clique);
	        return result;
	    }
	    if(clique.size()>=13) {
	        clique.addAll(p);
	        result.add(clique);
	        return result;
	      }
	    Set<Integer> intersectionPX = new HashSet<>(p);
	    intersectionPX.addAll(x);
	    int u = -1;
	    int max = 0;
	    for (Integer c : intersectionPX) {
	        int m = 0;
	        for (Integer v : adjacencyList.get(c)) {
	            if (intersectionPX.contains(v)) m++;
	        }
	        if (m > max){
	            max = m;
	            u = c;
	        }
	    }
	    Set<Integer> differenceP = new HashSet<>(p);
	    if (u != -1){
	        differenceP.removeAll(adjacencyList.get(u));
	    }
	    Iterator<Integer> it = differenceP.iterator();
	    while (it.hasNext()) {
	        Integer v = it.next();
	        Set<Integer> newCLique = new HashSet<>(clique);
	        newCLique.add(v);
	        Set<Integer> pIntersection = new HashSet<>(p);
	        Set<Integer> xIntersection = new HashSet<>(x);
	        pIntersection.retainAll(adjacencyList.get(v));
	        xIntersection.retainAll(adjacencyList.get(v));
	        result.addAll(findCliques(newCLique, pIntersection, xIntersection, adjacencyList));
	        p.remove(v);
	        x.add(v);
	    }
	    return result;
	}
	public static ArrayList<AminoAcid> Consensus(ArrayList<String> list){
		
		ArrayList<AminoAcid> FinalMotif  = new ArrayList<AminoAcid>();
		
		for(int i=0; i<7; i++) {
		Map<AminoAcid, Integer> inner1 = new HashMap<AminoAcid, Integer>();
		
		int vA=0;	int vR=0;	int vN=0;	int vD=0;	int vC=0;
		int vQ=0;	int vE=0;	int vG=0;	int vH=0;	int vI=0;
		int vL=0;	int vK=0;	int vM=0;	int vF=0;	int vP=0;
		int vS=0;	int vT=0;	int vW=0;	int vY=0;	int vV=0; int vDash=0;
		
		for (String s : list) {
			switch (s.charAt(i)) {
			case 'A':
				vA+=1;
				inner1.put(AminoAcid.A, vA);
				break;
			case 'R':
				vR+=1;
				inner1.put(AminoAcid.R, vR);
				break;
			case 'N':
				vN+=1;
				inner1.put(AminoAcid.N, vN);
				break;
			case 'D':
				vD+=1;
				inner1.put(AminoAcid.D, vD);
				break;
			case 'C':
				vC+=1;
				inner1.put(AminoAcid.C, vC);
				break;
			case 'Q':
				vQ+=1;
				inner1.put(AminoAcid.Q, vQ);
				break;
			case 'E':
				vE+=1;
				inner1.put(AminoAcid.E, vE);
				break;
			case 'G':
				vG+=1;
				inner1.put(AminoAcid.G, vG);
				break;
			case 'H':
				vH+=1;
				inner1.put(AminoAcid.H, vH);
				break;
			case 'I':
				vI+=1;
				inner1.put(AminoAcid.I, vI);
				break;
			case 'L':
				vL+=1;
				inner1.put(AminoAcid.L, vL);
				break;
			case 'K':
				vK+=1;
				inner1.put(AminoAcid.K, vK);
				break;
			case 'M':
				vM+=1;
				inner1.put(AminoAcid.M, vM);
				break;
			case 'F':
				vF+=1;
				inner1.put(AminoAcid.F, vF);
				break;
			case 'P':
				vP+=1;
				inner1.put(AminoAcid.P, vP);
				break;
			case 'S':
				vS+=1;
				inner1.put(AminoAcid.S, vS);
				break;
			case 'T':
				vT+=1;
				inner1.put(AminoAcid.T, vT);
				break;
			case 'W':
				vW+=1;
				inner1.put(AminoAcid.W, vW);
				break;
			case 'Y':
				vY+=1;
				inner1.put(AminoAcid.Y, vY);
				break;
			case 'V':
				vV+=1;
				inner1.put(AminoAcid.V, vV);
				break;
			case '-':
				vDash+=1;
				inner1.put(AminoAcid.Dash, vDash);
				break;
			}
		}
		 int maxValueInMap = (Collections.max(inner1.values()));  // This will return max value in the Hashmap
		   
		   for (Entry<AminoAcid, Integer> entry : inner1.entrySet()) {  // Itrate through hashmap
	         if (entry.getValue()==maxValueInMap) {
	            // System.out.println(entry.getKey());     // Print the key with max value
	             FinalMotif.add(entry.getKey());
	         }
	     }
	 }
		return FinalMotif;
	}
	
	public static int CountNum(int Cnt, ArrayList <Set<Integer>> Result) {
		int CntNew = 1;
		ArrayList<Integer> temp = new ArrayList<Integer>();
		for(int i= 0; i<Result.size(); i++) {
			temp.add(Result.get(i).size());
		}
		int max = Collections.max(temp);
		int countA = Collections.frequency(temp, max);
		Collections.sort(temp, Collections.reverseOrder());
		if (Cnt < countA) CntNew = countA;
		if (Cnt > countA) {
			temp.remove(max);
			max = Collections.max(temp);
			countA = Collections.frequency(temp, max);
			CntNew = Cnt + countA;
		}
		return CntNew;
	}
}
