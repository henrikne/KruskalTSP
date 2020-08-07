import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.Character.Subset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;


/**
 * KruskalTSP algorithm
 * Algorithm that uses modified Kruskal's algorithm to create open loop traveling salesman route
 * Also includes a version which uses randomization to get better result and basic Kruskal's algorithm implementation
 * @author Henrik Nenonen
 *
 */
public class KruskalTSP_Algorithm {
	
	
	/**
	 * Main
	 * @param args; filepath as string
	 */
	public static void main(String[] args) {
		
		//default file path
		String fileToRead = "distanceMatrix.txt";
		
		if (args.length != 0 && args != null)  {
			fileToRead = args[0];
			
		}
		//Reading input matrix
		ArrayList<ArrayList<Double>> inputMatrix = new ArrayList<ArrayList<Double>>();
		inputMatrix = readInputMatrix(fileToRead);
		System.out.println(inputMatrix);
		
		//Creating complete graph
		ArrayList<Node> nodes = new ArrayList<Node>();
		ArrayList<Edge> edgelist = new ArrayList<Edge>();
		Pair<ArrayList<Node>,ArrayList<Edge>> completegraph = CreateCompleteGraph(inputMatrix);
		nodes = completegraph.getValue1();
		edgelist = completegraph.getValue2();
		Collections.sort(edgelist);
		
		//Creating minimum spanning tree
		ArrayList<Edge>resultMST = CreateKruskalMST(nodes,edgelist);
		double resultlength = 0;
		
		//Printing results
		System.out.println("Printing MST");
		for (int i = 0; i < nodes.size()-1; i++) {
			resultlength = resultlength + resultMST.get(i).getLength();
			System.out.println("Nodes: "+resultMST.get(i).getStartPoint()+" ; "+ resultMST.get(i).getEndPoint());
			System.out.println(resultMST.get(i).getLength());
			
		}
		System.out.println(resultMST);
		System.out.println("Length: "+ resultlength);
		System.out.println();
		
		
		//########################################################################
		//Creating open loop traveling salesman route
		ArrayList<Edge>resultTSP = CreateKruskalTSP(nodes,edgelist);
		//ArrayList<Edge>resultTSP = CreateRandomizedKruskalTSP(graph);
		double resultTSPlength = 0;
		String resultTSProute = "";
		
		//Printing TSP edges
		System.out.println("Printing TSP route");
		for (int i = 0; i < nodes.size()-1; i++) {
			resultTSPlength = resultTSPlength + resultTSP.get(i).getLength();
			System.out.println("Nodes: "+resultTSP.get(i).getStartPoint()+" ; "+ resultTSP.get(i).getEndPoint());
			System.out.println(resultTSP.get(i).getLength());
			
		}
		//Finding leaf node from route
		int leafnode = 0;
		for (int j = 0; j < nodes.size(); j++) {
			if (nodes.get(j).getTimesintsp() == 1) {
				leafnode = j;
			}
		}
		System.out.println("leaf node: "+leafnode);
		
		//Constructing printable route from TSP edges
		ArrayList<Integer> resultroute = new ArrayList<Integer>();
		resultroute = constructRoute(resultTSP,leafnode);
		
		//Printing results
		System.out.println(resultTSP);
		System.out.println(resultTSProute);
		System.out.println("Result route: "+resultroute);
		System.out.println("TSP Length: "+ resultTSPlength);
		System.out.println("");
		
		//Reset how many times nodes are in TSP
		for (int k = 0; k < nodes.size(); k++) {
			nodes.get(k).resetTimesintsp();
		}
		
		//#######################################################################
		//Creating randomized TSP route
		int repeats = 100; //How many times we repeat making the randomized TSP route
		int rleafnode = 0;
		double bestRandomizedResultTSPlength = Double.MAX_VALUE;
		//String bestRandomizedResultTSProute = "";
		ArrayList<Edge>bestRandomizedResultTSP = new ArrayList<Edge>();
		
		for (int k = 0; k < repeats; k++) {
			
			//Reset nodes
			for (int j = 0; j < nodes.size(); j++) {
				nodes.get(j).resetTimesintsp();
			}
			
			ArrayList<Edge>randomizedResultTSP = CreateRandomizedKruskalTSP(nodes,edgelist);
			double randomizedResultTSPlength = 0;
			//String randomizedResultTSProute = "";
	
			//Printing TSP edges
			//System.out.println("Printing TSP route");
			for (int i = 0; i < nodes.size()-1; i++) {
				randomizedResultTSPlength = randomizedResultTSPlength + randomizedResultTSP.get(i).getLength();
				//System.out.println("Nodes: "+randomizedResultTSP.get(i).getStartPoint()+" ; "+ randomizedResultTSP.get(i).getEndPoint());
				//System.out.println(randomizedResultTSP.get(i).getLength());
	
			}
			if(randomizedResultTSPlength < bestRandomizedResultTSPlength) {
				bestRandomizedResultTSPlength = randomizedResultTSPlength;
				bestRandomizedResultTSP = randomizedResultTSP;
				
				//Finding leaf node from route
				for (int j = 0; j < nodes.size(); j++) {
					if (nodes.get(j).getTimesintsp() == 1) {
						rleafnode = j;
					}
				}
			}
			
			
		}
		
		System.out.println("randomized leaf node: "+rleafnode);

		//Constructing printable route from TSP edges
		ArrayList<Integer> randomizedResultroute = new ArrayList<Integer>();
		randomizedResultroute = constructRoute(bestRandomizedResultTSP,rleafnode);

		//Printing results
		System.out.println(bestRandomizedResultTSP);
		//System.out.println(randomizedResultTSProute);
		System.out.println("Randomized result route: "+randomizedResultroute);
		System.out.println("Randomized TSP Length: "+ bestRandomizedResultTSPlength);
		System.out.println();


	}
	
	
	public interface Pair<A,B>{
		public A getValue1();
		public B getValue2();
		
	}
	public static class Pair1<A,B> implements Pair<A,B>{

		private A value1;
		private B value2;
		
		public Pair1(A value1, B value2) {
			this.value1 = value1;
			this.value2 = value2;
		}
		
		@Override
		public A getValue1() {
			return value1;
		}

		@Override
		public B getValue2() {
			return value2;
		}
		
	}
	
	/** Reads file that has input matrix
	 * @param inputFilePath; location of the input file
	 * @return matrix with distances from every node to every other node
	 */
	public static ArrayList<ArrayList<Double>> readInputMatrix(String inputFilePath) {
		ArrayList<ArrayList<Double>> matrix=new ArrayList<ArrayList<Double>>();
		try {
			BufferedReader br=new BufferedReader(new FileReader(new File(inputFilePath)));
			String line;
			while((line=br.readLine())!=null){
				String[] lineItems = line.split(",");
				ArrayList<Double> row=new ArrayList<Double>();
				for(int j=0;j<lineItems.length;j++){
					row.add(Double.parseDouble(lineItems[j]));
				}
				matrix.add(row);
			}
			br.close();
		}
		catch (IOException e) {
			System.out.println("Can't read");
			e.printStackTrace();
			
			
		}
		//System.out.println(matrix);	
		return matrix;
	}
	
	/**
	 * Constructs printable route
	 * @param resultTSP; list of edges that are in TSP route
	 * @param currentnode; leaf node of TSP route
	 * @return route that has all nodes in route order
	 */
	private static ArrayList<Integer> constructRoute(ArrayList<Edge> resultTSP, int currentnode) {
		ArrayList<Integer> route = new ArrayList<Integer>();
		route.add(currentnode);
		
		while (route.size() < resultTSP.size()+1) {
			
			for(int i = 0; i < resultTSP.size(); i++) {
				
				if ((resultTSP.get(i).getStartPoint() == currentnode) && (route.contains(resultTSP.get(i).getEndPoint())==false)) {
					currentnode = resultTSP.get(i).getEndPoint();
					route.add(currentnode);		
				}
				
				else if ((resultTSP.get(i).getEndPoint() == currentnode) && (route.contains(resultTSP.get(i).getStartPoint())==false)) {
					currentnode = resultTSP.get(i).getStartPoint();
					route.add(currentnode);
					
				}	
			}	
		}
		return route;	
	}

	/**
	 * Creates complete graph from distance matrix
	 * @param inputMatrix; distance matrix (matrix that has distances from every node to every other node)
	 * @return 
	 */
	private static Pair<ArrayList<Node>, ArrayList<Edge>> CreateCompleteGraph(ArrayList<ArrayList<Double>> inputMatrix) {
		
		//edgelist.clear();
		ArrayList<Edge> edgelist = new ArrayList<Edge>();
		ArrayList<Node> graph = new ArrayList<Node>();
		ArrayList<String> uniqueEdges = new ArrayList<String>();
		for (int i = 0; i < inputMatrix.size(); i++ ) {
			Node nodenew = new Node();
			nodenew.setNumber(i);
			//nodenew.edges = null;
			System.out.println(nodenew.number);
			ArrayList<Double> distances = inputMatrix.get(i);
			for (int j = 0; j < distances.size(); j++) {
				Edge edgenew = new Edge();
				edgenew.setStartPoint(i);
				edgenew.setEndPoint(j);
				edgenew.setLength(distances.get(j));
				if (i != j) {
					nodenew.addEdge(edgenew);
				}
				if (!uniqueEdges.contains(edgenew.getStartPoint()+";"+edgenew.getEndPoint()) && !uniqueEdges.contains(edgenew.getEndPoint()+";"+ edgenew.getStartPoint()) && i != j) {					
					edgelist.add(edgenew);
					uniqueEdges.add(edgenew.getStartPoint()+";"+edgenew.getEndPoint());
					uniqueEdges.add(edgenew.getEndPoint()+";"+edgenew.getStartPoint());
				}
			}
		
			graph.add(nodenew);
			
		}
		
		Pair <ArrayList<Node>, ArrayList<Edge>> values = new Pair1 <ArrayList<Node>, ArrayList<Edge>>(graph, edgelist);
		return values; 
		
	}
	
	/**
	 * Creates minimum spanning tree using Kruskal's algorithm
	 * @param nodelist; 
	 * @return resulting MST as arraylist of edges
	 */
	private static ArrayList<Edge> CreateKruskalMST(ArrayList<Node> nodelist, ArrayList<Edge> edgelist) {
		
		ArrayList<Edge> result = new ArrayList<Edge>();
		System.out.println("edgelist: "+edgelist);
		
		System.out.println("sorted edgelist: "+edgelist);
		System.out.println("edgelist size: "+edgelist.size());
		
		//Create own set for each node and add them to list
		ArrayList<HashSet<Integer>> sets = new ArrayList<HashSet<Integer>>();
		System.out.println("Nodelist size: "+nodelist.size());
		int allnodessize = nodelist.size();
		for (int i = 0; i < allnodessize; i++) {
			HashSet<Integer> hash = new HashSet<Integer>();
			hash.add(nodelist.get(i).getNumber());
			sets.add(hash);
			
		}
		
		//While there is not enough edges in result loop through sorted edgelist
		int index = 0;
		//int sizeofedgelist = edgelist.size();
		while (result.size() < (nodelist.size()-1)){

			System.out.println("current subsets");
			System.out.println(sets);

			Edge testedge = edgelist.get(index);
			int startpoint = testedge.getStartPoint();
			int endpoint = testedge.getEndPoint();
			//System.out.println();
			//System.out.println("edgenodes: "+ startpoint + " ; "+endpoint);

			boolean foundinsameset = false;
			//Check if there is set where start and end are in same set
			for (int i = 0; i < sets.size(); i++) {
				if(sets.get(i).contains(startpoint) && sets.get(i).contains(endpoint)) {
					foundinsameset = true;
				}
			}

			//If nodes from tested edge are in same set, skip
			if(foundinsameset == true) {
				index++;
				System.out.println("edge skipped");
			}

			//Else mark nodes and add edge to result. Also remove edge from edgelist
			else {

				//Find proper indexes
				System.out.println("edge chosen");
				int indextoadd = 0;
				for(int i = 0; i < sets.size(); i++) {
					if (sets.get(i).contains(startpoint)) {
						indextoadd = i;
					}
				}
				int indexfrom = 0;
				for (int j = 0; j < sets.size(); j++) {
					if (sets.get(j).contains(endpoint)){
						indexfrom = j;
					}
				}
				sets.get(indextoadd).addAll(sets.get(indexfrom));
				sets.remove(indexfrom);
				//sets.get(endpoint).add(endpoint);
				result.add(testedge);
				index++;
				//edgelist.remove(index);
				//sizeofedgelist--;
				//break;
			}
		}
		return result;		
		
	}
	
	
	/**
	 * Creates TSP route using modified Kruskal's algorithm
	 * @param nodelist
	 * @return resulting arraylist of edges that contains TSP route 
	 */
	private static ArrayList<Edge> CreateKruskalTSP(ArrayList<Node> nodelist, ArrayList<Edge> edgelist) {
		ArrayList<Edge> result = new ArrayList<Edge>();
		System.out.println("edgelist: "+edgelist);
		
		//Collections.sort(edgelist);
		
		System.out.println("sorted edgelist: "+edgelist);
		System.out.println("edgelist size: "+edgelist.size());
		
		//Create own set for each node and add them to list
		ArrayList<HashSet<Integer>> sets = new ArrayList<HashSet<Integer>>();
		System.out.println("Nodelist size: "+nodelist.size());
		int allnodessize = nodelist.size();
		for (int i = 0; i < allnodessize; i++) {
			HashSet<Integer> hash = new HashSet<Integer>();
			hash.add(nodelist.get(i).getNumber());
			sets.add(hash);
			
		}
		
		
		//While there is not enough edges in result continue loop
		int index = 0;
		//int sizeofedgelist = edgelist.size();
		while (result.size() < (nodelist.size()-1)){

			System.out.println("current subsets");
			System.out.println(sets);

			//Choose shortest edge to test
			Edge testedge = edgelist.get(index);
			int startpoint = testedge.getStartPoint();
			int endpoint = testedge.getEndPoint();
			//System.out.println();
			System.out.println("edgenodes: "+ startpoint + " ; "+endpoint);

			boolean foundinsameset = false;
			//Check if there is set where start and end are in same set
			for (int i = 0; i < sets.size(); i++) {
				if(sets.get(i).contains(startpoint) && sets.get(i).contains(endpoint)) {
					foundinsameset = true;

				}
			}
			//Check if edge causes branching
			boolean branching = false;
			if (nodelist.get(startpoint).getTimesintsp() >= 2) {
				branching = true;
			}
			if (nodelist.get(endpoint).getTimesintsp() >= 2) {
				branching = true;
			}


			//If nodes from tested edge are in same set or edge would cause branching, skip
			if(foundinsameset == true || branching == true) {
				index++;
				System.out.println("edge skipped");
			}


			//Else mark nodes and add edge to result.
			else {

				//Find proper indexes
				System.out.println("edge chosen");
				int indextoadd = 0;
				for(int i = 0; i < sets.size(); i++) {
					if (sets.get(i).contains(startpoint)) {
						indextoadd = i;
					}
				}
				int indexfrom = 0;
				for (int j = 0; j < sets.size(); j++) {
					if (sets.get(j).contains(endpoint)){
						indexfrom = j;
					}
				}
				sets.get(indextoadd).addAll(sets.get(indexfrom));
				sets.remove(indexfrom);
				result.add(testedge);

				nodelist.get(startpoint).setTimesintsp();
				nodelist.get(endpoint).setTimesintsp();		
				index++;
			
			}
		}
		return result;
	}
	
	
	/**
	 * Creates TSP route using modified Kruskal's algorithm and randomization
	 * @param nodelist
	 * @return resulting arraylist of edges that contains TSP route 
	 */
	private static ArrayList<Edge> CreateRandomizedKruskalTSP(ArrayList<Node> nodelist, ArrayList<Edge> edgelist) {
		ArrayList<Edge> result = new ArrayList<Edge>();
		//System.out.println("edgelist: "+edgelist);
		//System.out.println("sorted edgelist: "+edgelist);
		//System.out.println("edgelist size: "+edgelist.size());
		
		//Create own set for each node and add them to list
		ArrayList<HashSet<Integer>> sets = new ArrayList<HashSet<Integer>>();
		//System.out.println("Nodelist size: "+nodelist.size());
		int allnodessize = nodelist.size();
		for (int i = 0; i < allnodessize; i++) {
			HashSet<Integer> hash = new HashSet<Integer>();
			hash.add(nodelist.get(i).getNumber());
			sets.add(hash);
			
		}
		
		ArrayList<Integer> availableEdges = new ArrayList<Integer>(Collections.nCopies(edgelist.size(), 0));
		//While there is not enough edges in result continue loop
		int index = 0;
		//int sizeofedgelist = edgelist.size();
		while (result.size() < (nodelist.size()-1)){

			//System.out.println("current subsets");
			//System.out.println(sets);

			//			System.out.println("result size: "+result.size());
			//			System.out.println("nodelist size: "+nodelist.size());
			
			//System.out.println(" edgelist size: "+sizeofedgelist);
			
			//Choose 5 shortest available edges and add their indexes to an arraylist
			ArrayList<Integer> chosenindexes = new ArrayList<Integer>();
			int testindex = 0;
			while (chosenindexes.size() < 5) {
				if (availableEdges.get(testindex) == 0) {
					chosenindexes.add(testindex);
				}
				testindex++;
			}
			//Choose randomly one of them
			int indexToTest = 0;
			Random rand = new Random();
			int randomnumber = rand.nextInt(5); 
			
			//If picking last edge to TSP route choose shortest one, because it is always the best pick
			if (result.size() == nodelist.size()-2) {
				indexToTest = chosenindexes.get(0);
			}
			else {
				indexToTest = chosenindexes.get(randomnumber);
			}
			
			Edge testedge = edgelist.get(indexToTest);
			
			int startpoint = testedge.getStartPoint();
			int endpoint = testedge.getEndPoint();
			//System.out.println();
			//System.out.println("edgenodes: "+ startpoint + " ; "+endpoint);

			//System.out.println("sets");
			//System.out.println(sets.get(startpoint));
			//System.out.println(sets.get(endpoint));
			boolean foundinsameset = false;
			
			//Check if there is set where start and end are in same set
			for (int i = 0; i < sets.size(); i++) {
				if(sets.get(i).contains(startpoint) && sets.get(i).contains(endpoint)) {
					foundinsameset = true;

				}
			}
			//Check if edge causes branching
			boolean branching = false;
			if (nodelist.get(startpoint).getTimesintsp() >= 2) {
				branching = true;
			}
			if (nodelist.get(endpoint).getTimesintsp() >= 2) {
				branching = true;
			}


			//If nodes from tested edge are in same set or edge would cause branching, skip
			if(foundinsameset == true || branching == true) {
				index++;
				availableEdges.set(indexToTest, 1);
				//System.out.println("edge skipped");
			}


			//Else mark nodes and add edge to result. 
			else {

				//Find proper indexes
				//System.out.println("edge chosen");
				availableEdges.set(indexToTest, 1);
				int indextoadd = 0;
				for(int i = 0; i < sets.size(); i++) {
					if (sets.get(i).contains(startpoint)) {
						indextoadd = i;
					}
				}
				int indexfrom = 0;
				for (int j = 0; j < sets.size(); j++) {
					if (sets.get(j).contains(endpoint)){
						indexfrom = j;
					}
				}
				sets.get(indextoadd).addAll(sets.get(indexfrom));
				sets.remove(indexfrom);
				//sets.get(endpoint).add(endpoint);
				result.add(testedge);
				//edgelist.remove(index);
				//sizeofedgelist--;


				nodelist.get(startpoint).setTimesintsp();
				nodelist.get(endpoint).setTimesintsp();		
				index++;
			
			}
		}
		return result;
	}
}


class Node implements Comparable<Object>{

	int number;
	int timesintsp = 0;
	ArrayList<Edge> edges = new ArrayList<Edge>();
	
	public int getNumber() {
		return number;
	}


	public void setNumber(int number) {
		this.number = number;
	}
	
	
	public int getTimesintsp() {
		return timesintsp;
	}


	public void setTimesintsp() {
		this.timesintsp++;
	}
	
	public void resetTimesintsp() {
		this.timesintsp = 0;
	}

	public void addEdge(Edge edge) {
		this.edges.add(edge);
	}
	//ArrayList<Double> coordinates;
	
	
	
	@Override
	public int compareTo(Object arg0) {
		// TODO Auto-generated method stub
		
		return 0;
	}
	
}

class Edge implements Comparable<Edge>{
	
	String id;
	int startpoint;
	int endpoint;
	double length;
	
	public String getId() {
		return id;
	}


	public void setId(String id) {
		this.id = id;
	}
	
	
	
	public int getStartPoint() {
		return startpoint;
	}


	public void setStartPoint(int startpoint) {
		this.startpoint = startpoint;
	}
	
	
	public int getEndPoint() {
		return endpoint;
	}


	public void setEndPoint(int endpoint) {
		this.endpoint = endpoint;
	}
	
	
	public double getLength() {
		return length;
	}


	public void setLength (double length) {
		this.length = length;
	}
	
	

	@Override
	public int compareTo(Edge o) {
		double compareLength = ((Edge) o).getLength();
		if (this.length > compareLength) {
			return 1;
		}
		else if (this.length < compareLength) {
			return -1;
		}
		else {
			return 0;
		}
	}
	

	
	
}


