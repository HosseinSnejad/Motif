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
import java.util.Set;
import java.util.Vector;

public class SequenceGeneration {
	public static enum AminoAcid {A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V, Dash};
	public static String aminoAcidLets = "ARNDCQEGHILKMFPSTWYV";
	
	public static void main(String [] args) throws FileNotFoundException {
		
		System.out.println();
		System.out.println("          "+ "Generating 10 Peptides of size 7:");
		System.out.println();
		
		String preDefSeq1 = "A--TG-Y";
		//String preDefSeq2 = "MR--E-Y";
	//	String preDefSeq3 = "Q-W--HN";
	//	String preDefSeq4 = "---LSTM";
	//	String preDefSeq5 = "-R-CI-V";
	//	String preDefSeq6 = "MK-F--P";
	//	String preDefSeq7 = "A-D-KT-";
		
		Vector sequences1 = sequenceGenerator(preDefSeq1, 1);
		//Vector sequences2 = sequenceGenerator(preDefSeq2, 5);
	//	Vector sequences3 = sequenceGenerator(preDefSeq3, 5);
	//	Vector sequences4 = sequenceGenerator(preDefSeq4, 5);
	//	Vector sequences5 = sequenceGenerator(preDefSeq5, 5);
	//	Vector sequences6 = sequenceGenerator(preDefSeq6, 5);
	//	Vector sequences7 = sequenceGenerator(preDefSeq7, 5);

		for(int i = 0;i < 1; i++)System.out.println("             " + sequences1.get(i));
		System.out.println();
		System.out.println();
		
		ArrayList<ArrayList<Integer>> indx = new ArrayList<ArrayList<Integer>>();
		ArrayList<Integer> t;
				
		Vector<String> SEQUENCES4 = new Vector<String>();
		
		ArrayList<ArrayList<Integer>> indxList = new ArrayList<ArrayList<Integer>>();

		for (int j = 15; j<= 120; j++){
			if (countSetBit(j) == 4)indxList.add(getSetIndxList(j));		
		}
		
		for(int i = 0; i < sequences1.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences1.get(i).toString(), indxList));
		}
	/*	for(int i = 0; i < sequences2.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences2.get(i).toString(), indxList));
		}
		for(int i = 0; i < sequences3.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences3.get(i).toString(), indxList));
		}*/
	/*	for(int i = 0; i < sequences4.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences4.get(i).toString(), indxList));
		}
		for(int i = 0; i < sequences5.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences5.get(i).toString(), indxList));
		}
		for(int i = 0; i < sequences6.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences6.get(i).toString(), indxList));
		}
		for(int i = 0; i < sequences7.size(); i++) {
			SEQUENCES4.addAll(kMerGenerator(sequences7.get(i).toString(), indxList));
		}
	*/	
		//  Writing to a File
		File file = new File("‎⁨Macintosh HD⁩ ▸ ⁨Users⁩ ▸ ⁨hosseinsaghaian⁩ ▸ ⁨eclipse-workspace⁩ ▸ ⁨1Motif⁩⁩");
		PrintStream out = new PrintStream(file);
		for(int i = 0; i < SEQUENCES4.size(); i++) {
			out.print(SEQUENCES4.get(i));
		}
		
		//Printing 4-Mers

		for(int m = 0; m < SEQUENCES4.size(); m++)System.out.println("             " + m + " : " + SEQUENCES4.get(m));
		System.out.println();
		
		
	}
	
	public static int countSetBit(int n) {
		int retVal = 0;
		while(n!=0) {
			retVal = retVal + (n & 1);
			n = n>>1;
		}
		return retVal;
	}

	public static ArrayList<Integer> getSetIndxList(Integer n){
		ArrayList<Integer> vect = new ArrayList<Integer> ();
		int temp = 0, cBit;
		
		while(n != 0) {
			cBit = n & 1;
			if (cBit == 1) {
				vect.add(temp);
			}
			else if (cBit==0 && temp == 0){
				vect.add(-7);}
			else {
				vect.add(-1 * temp);
			}
			n = n >> 1;
			temp++;
		}
		for(;temp < 7;temp++)vect.add(-1 * temp);
		return vect;
	}



	public static Vector sequenceGenerator(String preDefStr, int n ) {
		Vector retVal = new Vector();
		String temp;
		int randLetIndx;
		Random rand = new Random();
		boolean[] usedLet2 = new boolean[20];
		for(int i = 0; i < n; i++) {
			temp = "";
			
			for(int j = 0; j < preDefStr.length(); j++) {			
				
				if(preDefStr.charAt(j) == '-') {
					
					randLetIndx = rand.nextInt(20);
									
					temp += aminoAcidLets.charAt(randLetIndx);
									
				}
				else {
					temp = temp + preDefStr.charAt(j);
				}
			}
			retVal.add(temp);
		}
		return(retVal);
	}
	public static Vector kMerGenerator(String sequence, ArrayList<ArrayList<Integer>> indx) {
		Vector retVal = new Vector();
		String temp;
		for(int i=0; i<indx.size();i++) {
			temp = "";
			for(int k = 0; k < indx.get(i).size(); k++) {
				if(indx.get(i).get(k) >= 0) {
					temp += sequence.charAt(indx.get(i).get(k));
				}
				else {
					temp += '-';
				}
			}
			retVal.add(temp);
			
		}
		return(retVal);
	}

	public void indxGenerator() {
		ArrayList<int[]> indxList = new ArrayList<int[]>();
		for (int i = 1; i < 5; i++) {
			for (int j = 2; j > i && j < 6; j++) {
				for (int k = 3; k > j && k < 7; k++) {
					for (int l = 4; l > k; l++) {
						int [] array = {i, j, k, l};
						indxList.add(array);
					}
				}
			}
		}
		System.out.print("{");
		for (int p = 0; p< indxList.size(); p++) {
			System.out.print(Arrays.toString(indxList.get(p)));
		}
		System.out.print("}");
	}

	public static ArrayList<Integer> getNullIndxList(ArrayList<Integer> List){
		ArrayList<Integer> NullIndx = new ArrayList<Integer>();
		for (int i = 1; i<=3;i++) {
			switch (List.get(i)-List.get(i-1)){
			case 1:
				break;
		case 2:
			NullIndx.add(List.get(i-1)+1);
			break;
		case 3:
			NullIndx.add(List.get(i-1)+1);
			NullIndx.add(List.get(i-1)+2);
			break;
		case 4:
			NullIndx.add(List.get(i-1)+1);
			NullIndx.add(List.get(i-1)+2);
			NullIndx.add(List.get(i-1)+3);
			break;
			}
		}
		return NullIndx;
	}

	public static ArrayList<Integer> getSetIndxList1(Integer n){
		ArrayList<Integer> vect = new ArrayList<Integer> ();
		int temp = 0, cBit;
		while(n != 0) {
			cBit = n & 1;
			if (cBit == 1) {
				vect.add(temp);
			}
			n = n >> 1;
			temp++;
		}
		return vect;
	}

	public static int[][] AdjacencyMatrix(Vector<String> seq){
		int[][] adjMatrix = new int [350][350];
		for(int i=0; i<seq.size(); i++) {
			for(int j=i+1; j<seq.size()-1; j++) {
				int score=0;
				for(int k = 0; k < 7; k++) {
					if(seq.get(i).charAt(k) == seq.get(j).charAt(k)) {
						score+=1;
					}
					
			}
				if (score >= 4)adjMatrix[i][j] = 1;
				}
			}
		return adjMatrix;
		}

	public static ArrayList<Set<Integer>> AdjList(Vector<String> seq){
		ArrayList<Set<Integer>> adjacencyList = new ArrayList();
		for(int i=0; i<seq.size(); i++) {
			adjacencyList.add(new HashSet());
		}
		for(int i=0; i<seq.size(); i++) 
		{
			for(int j = 0; j < seq.size(); j++) 
			{
				if(i == j)continue;
				int score=0;
				for(int k = 0; k < 7; k++) 
				{
					if( seq.get(i).charAt(k) == seq.get(j).charAt(k) && seq.get(i).charAt(k) !=0 ) 
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
	             System.out.println(entry.getKey());     // Print the key with max value
	             FinalMotif.add(entry.getKey());
	         }
	     }
	 }
		return FinalMotif;
	}
}
