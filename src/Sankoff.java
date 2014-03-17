 /**
 * Created by Gordon on 2/11/14.
 */
import java.io.*;
import java.util.*;

 /* outline
  *Open up Stockholm and get data u
  *open up tree and get data about phylogenetic tree --> put into actual tree format
  *bottom down
  *top up
  */
public class Sankoff {
    public static ArrayList<PhyloTreeNode> sequenceVariations; //ArrayList PhyloTreeNodes that represent a gene variation in the family
    public static String ss_Cons = null; //consensus secondary structure
    public static String rf = null; //sequence alignment???
    public static PhyloTree tree; //phylogenetic tree of the gene family
    public static String name; //name of the family we are working on
    public static String alphabet = "ABCDEFGHIJKLMNOPQRSTUV"; //alphabet so that each new node gets a 1 letter unique name
    public static int alphCounter = 0; //counter so that each new node gets a unique name
    private static int[][] cost;
    private static final int INF = 10000;

    public static void stockholmParse(String stockholmFile) throws IOException{
        sequenceVariations = new ArrayList<PhyloTreeNode>();
        String line = null;
        int sequenceCounter = 0; //counter to align the sequences to the PhyloTreeNodes with the correct names
        name = ""; //name of the family we are working on
        try{
            FileReader fileReader = new FileReader(stockholmFile);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            while((line = bufferedReader.readLine()) != null){
                if(line.contains("#=GF DE")){
                    String[] parts = line.split("\\s{1,}");
                    for(int i = 2; i < parts.length; i++){
                        name = name.concat(parts[i]); //line with the
                        name = name.concat(" ");
                    }
                    name = name.substring(0, name.length() - 1);
                }
                if(line.contains("#=GC SS_cons")){
                    String[] parts = line.split("\\s{1,}");
                    ss_Cons = parts[2];
                }
                if(line.contains("#=GC RF")){
                    String[] parts = line.split("\\s{1,}");
                    rf = parts[2];
                }
                if(line.contains("#=GS")){
                    String[] parts = line.split("\\s{1,}");
                    PhyloTreeNode current = new PhyloTreeNode(parts[3], null);
                    sequenceVariations.add(current);
                }
                if(!line.startsWith("#") && !line.equals("") && sequenceCounter < sequenceVariations.size()){
                    String[] parts = line.split("\\s{1,}");
                    sequenceVariations.get(sequenceCounter).setSequence(parts[1]);
                    sequenceCounter ++;
                }
            }
        }
        catch(FileNotFoundException ex) {
            System.out.println("no file homie");
        }
        catch(IOException ex){
            System.out.println("error reading the file");
        }
        //System.out.println(name);
        //System.out.println(ss_Cons);
        //System.out.println(rf);
        //for(int i = 0; i < sequenceVariations.size(); i++){
        //System.out.println(sequenceVariations.get(i).getName());
        //System.out.println(sequenceVariations.get(i).getSequence());
        //}
    }

    //identify the most inner bracket --> find the first )
    //find the ( that goes with it
    //combine the 2 of them to create a single node, the root, labeled A, will be their parents
    //recursively move to the next layer, where you combine the the current tree and the new sequence, by creating a new root, with its
    //children being the old root, and the new sequence
    //if there are 2 genes in a single layer at any point,
    //create a parent for them, treat the parent node like a single node
    //testd
    public static void phyloTreeCreator(String treeFile) throws IOException{
        String line = null;
        int curOpen;
        int curClosed;
        tree = new PhyloTree(alphabet.substring(alphCounter, alphCounter + 1), "tbd"); //creates root with name A
        alphCounter++;
        tree.setConsensusSequence(ss_Cons);
        try{
            FileReader fileReader = new FileReader(treeFile);
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            while((line = bufferedReader.readLine()) != null){

                //Next step is to figure out how many genes are in-between
                //find out which genes are in between
                //if 2, add them together and treat parent as if there was a single gene
                //create a new root, make the previous root and the gene the children of the new root
                while(line.contains("(") || line.contains(")")){ //after dealing with a layer, swap out () so when done there are none left
                    curClosed = line.indexOf(")"); //finds first instance of ) in the string
                    curOpen = findClosestOpen(curClosed, line); //method returns closest ( bracket in the string
                    String layer = line.substring(curOpen, curClosed + 1); //creates layer, which is the most interior bracket set
                    int numberofGenesInLayer = numberOfSequences(layer); //finds if there is 1 or 2 genes in the layer
                    PhyloTreeNode newNode = null;

                    //if there are 2 nodes in the layer, we are going to combine them and create a parent node
                    //if the root has < 2 children, just attach this new node
                    //if the root has 2 children, then we need to create a new root node which will have children
                    //that are our combined nodes and the old root
                    if(numberofGenesInLayer == 2){
                        newNode = twoLayer(layer);
                        if(tree.getRoot().getChildren().size() == 2){
                            PhyloTreeNode newRoot = new PhyloTreeNode(alphabet.substring(alphCounter, alphCounter + 1), "tbd");
                            newRoot.addChild(newNode);
                            tree.addNode(newRoot);
                            alphCounter++;
                        }
                        else tree.getRoot().addChild(newNode);
                    }
                    //if the layer only has 1 gene
                    //if the root has < 2 children, just add the node from the layer as a child to the root
                    //if the root has 2 children, create a new root with children that are our new node, and the old root
                    if(numberofGenesInLayer == 1){
                        String curName;
                        for(int i = 0; i < sequenceVariations.size(); i++){
                            curName = sequenceVariations.get(i).getName();
                            if(layer.indexOf(curName) > 0){
                                newNode = sequenceVariations.get(i);
                            }
                        }
                        if(tree.getRoot().getChildren().size() == 2){
                            PhyloTreeNode newRoot = new PhyloTreeNode(alphabet.substring(alphCounter, alphCounter + 1), "tbd");
                            newRoot.addChild(newNode);
                            tree.addNode(newRoot);
                            alphCounter++;
                        }
                        else tree.getRoot().addChild(newNode);
                    }
                    //method to cut out the old layer
                    if(curOpen != 0){
                        String lineFirstHalf = line.substring(0, curOpen);
                        String lineSecondHalf = line.substring(curClosed + 1);
                        line = lineFirstHalf.concat(lineSecondHalf);
                    }
                    else break;
                }
            }
        }
        catch(FileNotFoundException ex) {
            System.out.println("no file homie1");
        }
        catch(IOException ex){
            System.out.println("error reading the file");
        }
    }
    //method that takes a layer of 2 and combines the 2 nodes into a new parent node
    private static PhyloTreeNode twoLayer(String layer) {
        String curName;
        PhyloTreeNode parent = new PhyloTreeNode(alphabet.substring(alphCounter, alphCounter + 1), "tbd");
        alphCounter++;
        for(int i = 0; i < sequenceVariations.size(); i++){
            curName = sequenceVariations.get(i).getName();
            if(layer.indexOf(curName) > 0){
                parent.addChild(sequenceVariations.get(i));
            }
        }
        return parent;
    }
    //figures out if the layer is a one or 2 layer
    private static int numberOfSequences(String layer){
        int counter = 0;
        String curName;
        for(int i = 0; i < sequenceVariations.size(); i++){
            curName = sequenceVariations.get(i).getName();
            if(layer.indexOf(curName) > 0) counter++;
        }
        return counter;
    }

    //Finds the closest open bracket to our current closed bracker
    private static int findClosestOpen(int closedLocation, String line){
        int curClosestOpen = 0;
        int nextOpen = 0;
        while(nextOpen < closedLocation && nextOpen != -1){
            curClosestOpen = nextOpen;
            nextOpen = line.indexOf("(", curClosestOpen + 1);
        }
        return curClosestOpen;
    }
    //BFS print out of the tree
    private static void printTree(PhyloTree tree){
        PhyloTreeNode current = null;
        Queue<PhyloTreeNode> s = new LinkedList<PhyloTreeNode>();
        s.add(tree.getRoot());
        while(!s.isEmpty()){
            current = s.poll();
            System.out.print(current.getName());
            if(current.getParent() != null) System.out.print("     " + current.getParent().getName());
            else System.out.print("       no parent");
            System.out.println("     " + current.getSequence());
            if(current.getChildren() != null){
                for(int i = 0; i < current.getChildren().size(); i++){
                    s.add(current.getChildren().get(i));
                }
            }
        }
    }


    private static void sequenceMod(){
        for(int i = 0; i < sequenceVariations.size(); i++){
            String seq = sequenceVariations.get(i).getSequence();
            seq = seq.replace("", ",");
            seq = seq.substring(1, seq.lastIndexOf(","));
            sequenceVariations.get(i).setSequence(seq);
        }
    }

    private static int bottomUp(PhyloTreeNode n) {
        PhyloTreeNode cur = null;
        int parsimonyScore = 0;
        for(int i = 0; i < n.getChildren().size(); i++){
            cur = n.getChildren().get(i);
            if(cur.getSequence().equals("tbd")){
                parsimonyScore += bottomUp(cur);
            }
        }
        String newSeq = "";
        if(n.getChildren().size() == 1) n.setSequence(n.getChildren().get(0).getSequence());
        if(n.getChildren().size() == 2){
            String s1 = n.getChildren().get(0).getSequence();
            String s2 = n.getChildren().get(1).getSequence();
            StringTokenizer t1 = new StringTokenizer(s1, ",");
            StringTokenizer t2 = new StringTokenizer(s2, ",");
            int tokensNum = t1.countTokens();
            for(int j = 0; j < tokensNum; j++){
                String c1 = t1.nextToken();
                String c2 = t2.nextToken();
                if(c1.equals(c2)){
                    newSeq = newSeq.concat(c1);
                    newSeq = newSeq.concat(",");
                }
                else{
                    boolean oneMatch = false;
                    char[] singleBase1 = c1.toCharArray();
                    char[] singleBase2 = c2.toCharArray();
                    for(int i = 0; i < singleBase1.length; i++){
                        for(int k = 0; k < singleBase2.length; k++){
                            if (singleBase1[i] == singleBase2[k]){
                                oneMatch = true;
                                newSeq = newSeq.concat(java.lang.Character.toString(singleBase1[i]));
                            }
                        }
                    }
                    if(!oneMatch){
                        for(int i = 0; i < singleBase1.length; i++){
                            newSeq = newSeq.concat(java.lang.Character.toString(singleBase1[i]));
                        }
                        for(int k = 0; k < singleBase2.length; k++){
                            newSeq = newSeq.concat(java.lang.Character.toString(singleBase2[k]));
                        }
                        parsimonyScore++;
                    }
                    newSeq = newSeq.concat(",");
                }
            }
        }
        n.setSequence(newSeq);
        return parsimonyScore;
    }
    //just recursivley go through all the nodes
    //go through the sequences and if there are ever any more than one letter, refer to the parent
    //if the parent is one of the options, pick that one, if that just pick the first one
    private static void topDown(PhyloTreeNode node){
        if(node.getChildren().size() != 0){
            for(int i = 0; i < node.getChildren().size(); i++){
                topDown(node.getChildren().get(i));
            }
        }
        boolean parent = false;
        StringTokenizer parentTokens = null;
        if(node.getParent() != null){
            parent = true;
            parentTokens = new StringTokenizer(node.getParent().getSequence(), ",");
        }
        StringTokenizer nodeTokens = new StringTokenizer(node.getSequence(), ",");
        String newSeq = "";
        int tokenCount = nodeTokens.countTokens();
        String curParentToken = null;
        for(int i = 0; i < tokenCount; i++){
            String curNodeToken = nodeTokens.nextToken();
            if(parent) curParentToken = parentTokens.nextToken();
            if(curNodeToken.length() > 1){
                if(parent){
                    boolean match = false;
                    for(int j = 0; j < curNodeToken.length(); j++){
                        if(curNodeToken.charAt(j) == curParentToken.charAt(0)){
                            newSeq = newSeq.concat(curParentToken);
                            match = true;
                            break;
                        }
                    }
                    if(!match){
                        newSeq = newSeq.concat(curNodeToken.substring(0,1));
                    }
                }
                else newSeq = newSeq.concat(curNodeToken.substring(0,1));
            }
            else newSeq = newSeq.concat(curNodeToken.substring(0,1));


        }
        node.setSequence(newSeq);
    }
    public static int sankoff(PhyloTree tree){
        int min = -INF;
        int cur = 0;
        int totalScore = 0;
        Collection<String> bases = new ArrayList<String>();
        bases.add("A");
        bases.add("C");
        bases.add("G");
        bases.add("U");
        bases.add(".");
        String bestBase = "";
        String curBase = "";
        PhyloTreeNode curNode = tree.getRoot();
        while(curNode.getChildren().size() == 2) curNode = curNode.getChildren().get(1);
        String seq = curNode.getSequence();
        seq = seq.replace(",","");
        int seqLength = seq.length();
        for(int i = 0; i < seqLength; i++){
            Iterator<String> it = bases.iterator();
            while(it.hasNext()){
                curBase = it.next();
                cur = sankoffRecursion(tree.getRoot(), curBase, i);
                if(cur > min){
                    min = cur;
                    bestBase = curBase;
                }
            }
            if(i == 0) tree.getRoot().setSequence(bestBase);
            else tree.getRoot().setSequence(tree.getRoot().getSequence().concat(bestBase));
            totalScore += min;
            min = -INF;
        }
        System.out.println(totalScore);
        Queue<PhyloTreeNode> q = new LinkedList<PhyloTreeNode>();
        q.add(tree.getRoot());
        PhyloTreeNode curN;
        while(! q.isEmpty()){
            curN = q.poll();
            for(int i = 0; i < curN.getChildren().size(); i++){
                q.add(curN.getChildren().get(i));
            }
            String curSeq = curN.getSequence();
            curSeq = curSeq.replace(",","");
            curSeq = curSeq.replace("", ",");
            if(curSeq.startsWith(",")) curSeq = curSeq.replaceFirst(",", "");
            if(curSeq.endsWith(",")) curSeq = curSeq.substring(0,curSeq.length() -1);
            curN.setSequence(curSeq);
        }

        return totalScore;
    }

    private static int sankoffRecursion(PhyloTreeNode node, String base, int pos){
        int totalMax = -INF;
        int totalSum = 0;
        String bestBase = " ";
        Collection<String> bases = new ArrayList<String>();
        bases.add("A");
        bases.add("C");
        bases.add("G");
        bases.add("U");
        bases.add(".");
        //lets pretend that the base in the current position is the base passed
        Iterator<String> itL = bases.iterator();
        Iterator<String> itR = bases.iterator();
        String curL;
        String curR;
        String seqMod = node.getSequence().replace(",", "");
        String[] seqArray = seqMod.split("");
        String bestL = "";
        String bestR = "";
        //if we are at a leaf
        //if our pretend base matches the actual base that is there
        //score of 0, ow score of -INF
        if(node.getChildren().size() < 2){
           if(seqArray[pos+1].equals(base)){
               return 0;
           }
           else{
               return -INF;
           }
        }    //the score to go from our current pretend base to a new base
        //plus the sankoff of running the child
        else{
            int maxL = -INF * 2;
            int sum;
            while(itL.hasNext()){
                curL = itL.next();
                sum = cost(base, curL);
                sum += sankoffRecursion(node.getChildren().get(0), curL,pos);
                if(sum > maxL){
                    maxL = sum;
                    bestL = curL;
                }
            }
            int maxR = -INF * 2;
            while(itR.hasNext()){
                curR = itR.next();
                sum = cost(base, curR);
                sum += sankoffRecursion(node.getChildren().get(1), curR, pos);
                if(sum > maxR){
                    maxR = sum;
                    bestR = curR;
                }
            }
            totalMax = maxR + maxL;

            //non-essential code to get sequences at each letter
            if(pos == 0 && node.getChildren().get(0).getChildren().size() > 1){
                node.getChildren().get(0).setSequence(bestL);
                node.getChildren().get(1).setSequence(bestR);
            }
            if(pos > 0 && node.getChildren().get(0).getChildren().size() > 1){
                PhyloTreeNode leftChild = node.getChildren().get(0);
                PhyloTreeNode rightChild = node.getChildren().get(1);
                if(!(pos < leftChild.getSequence().length())){
                    leftChild.setSequence(leftChild.getSequence().concat(bestL));
                    rightChild.setSequence(rightChild.getSequence().concat(bestR));
                }
                else{
                    String[] leftSeq = leftChild.getSequence().split("");
                    String[] rightSeq = rightChild.getSequence().split("");
                    leftSeq[pos + 1] = bestL;
                    rightSeq[pos + 1] = bestR;
                    String newL = "";
                    for(String str: leftSeq)newL += str;
                    String newR = "";
                    for(String str: rightSeq)newR += str;
                    leftChild.setSequence(newL);
                    rightChild.setSequence(newR);
                }
            }
            //end of non-essential code
            PhyloTreeNode current = node;
            while(current.getChildren().size() > 1){
                current = current.getChildren().get(0);
            }
            int seqLength = current.getSequence().length();
            ArrayList<Integer> differences = new ArrayList<Integer>();
            if(node.getPontSequence().length() == seqLength){
                String oldSeq = node.getSequence();
                String pontSeq = node.getPontSequence();
                for(int i = 0; i < seqLength; i++){
                    if(oldSeq.charAt(i) == pontSeq.charAt(i)) differences.add(i);
                }

            }

        }
        return totalMax;
    }
        //if we aren't at a leaf, for our pretend base
        //the score for our base is dependent on 2 criteria

    private static int cost(String b1, String b2){
        int b1Int = 0;
        int b2Int = 0;
        switch (b1.charAt(0)){
            case 'A':
                b1Int = 0;
                break;
            case 'C':
                b1Int = 1;
                break;
            case 'G':
                b1Int = 2;
                break;
            case 'U':
                b1Int = 3;
                break;
            case '.':
                b1Int = 4;
                break;
        }
        switch (b2.charAt(0)){
            case 'A':
                b2Int = 0;
                break;
            case 'C':
                b2Int = 1;
                break;
            case 'G':
                b2Int = 2;
                break;
            case 'U':
                b2Int = 3;
                break;
            case '.':
                b2Int = 4;
                break;
        }
        return cost[b1Int][b2Int];
    }

    public static void rnaFold(PhyloTree tree) throws IOException {
        String command = "C:\\Users\\Gordon\\Dropbox\\Winter2014\\Comp401\\ViennaRNAPackage\\rnaFold.exe";
        BufferedReader inp;
        BufferedWriter out;
        ProcessBuilder builder = new ProcessBuilder(command);
        builder.redirectErrorStream(true);
        Process p = builder.start();
        InputStream ips = p.getInputStream();
        OutputStream ops = p.getOutputStream();
        inp = new BufferedReader(new InputStreamReader(ips));
        out = new BufferedWriter(new OutputStreamWriter(ops));
        Queue<PhyloTreeNode> q = new LinkedList<PhyloTreeNode>();
        q.add(tree.getRoot());
        PhyloTreeNode cur;
        while(! q.isEmpty()){
            cur = q.poll();
            for(int i = 0; i < cur.getChildren().size(); i++){
                q.add(cur.getChildren().get(i));
            }
            String seq = cur.getSequence();
            seq = seq.replace(".", "");
            seq = seq.concat("\n");
            out.write(seq);
            out.flush();
            String line;
            int i = 0;
            while(i < 2 && ( line = inp.readLine()) != null){
                if(i == 1){
                    int lastClosed = line.lastIndexOf(")");
                    int lastOpen = findClosestOpen(lastClosed, line);
                    String energy = line.substring(lastOpen);
                    line = line.substring(0, lastOpen);
                    cur.setFolding(line);
                    energy = energy.replace("(", "");
                    energy = energy.replace(")", "");
                    energy =  energy.trim();
                    cur.setEnergy(Double.parseDouble(energy));
                }
                i++;
            }
            System.out.println(cur.getFolding() + "   " + Double.toString(cur.getEnergy()));
        }
        p.destroy();

    }
    public static void calcDistancesFromConsensus(PhyloTree tree) throws IOException {
        String command = "C:\\Users\\Gordon\\Dropbox\\Winter2014\\Comp401\\ViennaRNAPackage\\rnaDistance.exe";
        BufferedReader inp;
        BufferedWriter out;
        ProcessBuilder builder = new ProcessBuilder(command);
        builder.redirectErrorStream(true);
        Process p = builder.start();
        InputStream ips = p.getInputStream();
        OutputStream ops = p.getOutputStream();
        inp = new BufferedReader(new InputStreamReader(ips));
        out = new BufferedWriter(new OutputStreamWriter(ops));
        Queue<PhyloTreeNode> q = new LinkedList<PhyloTreeNode>();
        q.add(tree.getRoot());
        PhyloTreeNode cur;
        while(! q.isEmpty()){
            cur = q.poll();
            for(int i = 0; i < cur.getChildren().size(); i++){
                q.add(cur.getChildren().get(i));
            }
            String con = tree.getConsensusSequence().concat("\n");
            con = con.replace("<", "(");
            con = con.replace(">", ")");
            String curFolding = cur.getFolding().concat("\n");
            out.write(con);
            out.flush();
            out.write(curFolding);
            out.flush();
            String line;
            int i = 0;
            while(i < 1 && ( line = inp.readLine()) != null){
                //System.out.println(line);
                line = line.substring(2);
                line = line.trim();
                cur.setDistanceFromConsensus(Integer.parseInt(line));
                i++;
            }
            System.out.println(cur.getDistanceFromConsensus());
        }


    }
     public static int keepStructureRecursive(PhyloTreeNode node, String base, int pos){
         int totalMax = -INF;
         int totalSum = 0;
         String bestBase = " ";
         Collection<String> bases = new ArrayList<String>();
         bases.add("A");
         bases.add("C");
         bases.add("G");
         bases.add("U");
         bases.add(".");
         //lets pretend that the base in the current position is the base passed
         Iterator<String> itL = bases.iterator();
         Iterator<String> itR = bases.iterator();
         String curL;
         String curR;
         String seqMod = node.getSequence().replace(",", "");
         String[] seqArray = seqMod.split("");
         String bestL = "";
         String bestR = "";
         //if we are at a leaf
         //if our pretend base matches the actual base that is there
         //score of 0, ow score of -INF
         if(node.getChildren().size() < 2){
             if(seqArray[pos+1].equals(base)){
                 return 0;
             }
             else{
                 return -INF;
             }
         }    //the score to go from our current pretend base to a new base
         //plus the sankoff of running the child
         else{
             int maxL = -INF * 2;
             int sum;
             while(itL.hasNext()){
                 curL = itL.next();
                 sum = cost(base, curL);
                 sum += keepStructureRecursive(node.getChildren().get(0), curL, pos);
                 if(sum > maxL){
                     maxL = sum;
                     bestL = curL;
                 }
             }
             int maxR = -INF * 2;
             while(itR.hasNext()){
                 curR = itR.next();
                 sum = cost(base, curR);
                 sum += keepStructureRecursive(node.getChildren().get(1), curR, pos);
                 if(sum > maxR){
                     maxR = sum;
                     bestR = curR;
                 }
             }
             totalMax = maxR + maxL;

             //non-essential code to get sequences at each letter
             if(pos == 0 && node.getChildren().get(0).getChildren().size() > 1){
                 node.getChildren().get(0).setPontSequence(bestL);
                 node.getChildren().get(1).setPontSequence(bestR);
             }
             if(pos > 0 && node.getChildren().get(0).getChildren().size() > 1){
                 PhyloTreeNode leftChild = node.getChildren().get(0);
                 PhyloTreeNode rightChild = node.getChildren().get(1);
                 if(!(pos < leftChild.getPontSequence().length())){
                     leftChild.setPontSequence(leftChild.getPontSequence().concat(bestL));
                     rightChild.setPontSequence(rightChild.getPontSequence().concat(bestR));
                 }
                 else{
                     String[] leftSeq = leftChild.getPontSequence().split("");
                     String[] rightSeq = rightChild.getPontSequence().split("");
                     leftSeq[pos + 1] = bestL;
                     rightSeq[pos + 1] = bestR;
                     String newL = "";
                     for(String str: leftSeq)newL += str;
                     String newR = "";
                     for(String str: rightSeq)newR += str;
                     leftChild.setPontSequence(newL);
                     rightChild.setPontSequence(newR);
                 }
             }

         }
         return totalMax;
     }

     public static void keepSequence(PhyloTree tree){
         PhyloTreeNode curNode;
         int seqLength;
         int totalScore = 0;
         Collection<String> bases = new ArrayList<String>();
         bases.add("A");
         bases.add("C");
         bases.add("G");
         bases.add("U");
         bases.add(".");
         curNode = tree.getRoot();
         while((curNode.getChildren().size()) > 0 ) curNode = curNode.getChildren().get(0);
         String seq = curNode.getSequence();
         seq = seq.replace(",","");
         seqLength = seq.length();
         for(int i = 0; i < seqLength; i++){
             Iterator<String> it = bases.iterator();
             int bestScore = -INF;
             String bestBase = "";
             while(it.hasNext()){
                 String curBase = it.next();
                 int curBaseScore = keepStructureRecursive(tree.getRoot(), curBase, i);
                 if(curBaseScore > bestScore){
                    bestScore = curBaseScore;
                    bestBase = curBase;
                 }
             }
             tree.getRoot().setPontSequence(tree.getRoot().getPontSequence().concat(bestBase));
             totalScore += bestScore;
         }
         Queue<PhyloTreeNode> q = new LinkedList<PhyloTreeNode>();
         q.add(tree.getRoot());
         PhyloTreeNode curN;
         while(! q.isEmpty()){
             curN = q.poll();
             for(int i = 0; i < curN.getChildren().size(); i++){
                 q.add(curN.getChildren().get(i));
             }
             String curSeq = curN.getPontSequence();
             curSeq = curSeq.replace(",","");
             curSeq = curSeq.replace("", ",");
             if(curSeq.startsWith(",")) curSeq = curSeq.replaceFirst(",", "");
             if(curSeq.endsWith(",")) curSeq = curSeq.substring(0,curSeq.length() -1);
             curN.setPontSequence(curSeq);
         }


     }



    public static void main(String Args[]){
        try {
            stockholmParse(Args[0]);
            phyloTreeCreator(Args[1]);
            //printTree(tree);
            sequenceMod();
            //System.out.println(bottomUp(tree.getRoot()));
            //printTree(tree);
            //topDown(tree.getRoot());
            //printTree(tree);
            //rnaFold(tree);
            //calcDistancesFromConsensus(tree);
            cost = new int[5][5];
            for(int i = 0; i < 5; i++){
                if(i != 4){
                    cost[i][4] = -2;
                    cost[4][i] = -2;
                }
                cost[i][i] = 0;
            }
            //get across the diagonal correct
            cost[0][1] = -2;
            cost[0][2] = -1;
            cost[0][3] = -2;
            cost[1][2] = -2;
            cost[1][3] = -1;
            cost[2][3] = -2;
            cost[1][0] = -2;
            cost[2][0] = -1;
            cost[3][0] = -2;
            cost[2][1] = -2;
            cost[3][1] = -1;
            cost[3][2] = -2;
            int j = 0;
            sankoff(tree);
            printTree(tree);
            keepSequence(tree);
            Queue<PhyloTreeNode> q = new LinkedList<PhyloTreeNode>();
            PhyloTreeNode curNode = tree.getRoot();
            q.add(curNode);
            while(!q.isEmpty()){
                curNode = q.poll();
                for(int k = 0; k < curNode.getChildren().size(); k++){
                    q.add(curNode.getChildren().get(k));
                }
                System.out.println(curNode.getPontSequence());
            }







        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}