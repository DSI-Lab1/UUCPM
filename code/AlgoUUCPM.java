import java.io.*;
import java.util.*;
import java.util.Map.Entry;
/**
 * Copyright (C), 2022-2023, Jinan University
 * FileName: AlgoUUCPM
 * Author:   Zefeng Chen
 * Date:     2020/1/25 10:20
 * Description: FUCPM algorithm.
 * Pruning strategy: Recurrent PIWU + PIEU.
 * Data structure: Sequence Infomation List and Instance List.
 */
public class AlgoUUCPM {

    /** the time the algorithm started */
    double startTimestamp = 0;
    /** the time the algorithm terminated */
    long endTimestamp = 0;
    /** the number of patterns generated */
    int patternCount = 0;

    /** writer to write the output file **/
    BufferedWriter writer = null;

    /** buffer for storing the current pattern that is mined when performing mining
     * the idea is to always reuse the same buffer to reduce memory usage. **/
    final int BUFFERS_SIZE = 2000;
    private int[] patternBuffer = null;

    /** if true, debugging information will be shown in the console */
    //1:PIWU, 2:recurrent PIWU, 3:SIL, 4:InstanceChain of 1-seq, 5:utility of 1-seq, 6:projected InstanceChain, 7:prefix utility, 8:writeout function
    final int DEBUG = 0;

    /** if true, save result to file in a format that is easier to read by humans **/
    final boolean SAVE_RESULT_EASIER_TO_READ_FORMAT = false;

    /** the minUtility threshold **/
    double minUtility = 0;

    /** max pattern length **/
    int maxPatternLength = 1000;

    /** the input file path **/
    String input;
    
    boolean GUIRPstrategy = true;
    boolean UEIPstrategy = false;
    
    boolean item_Level = false; 
    //if it is true, the UUCPM supports database with probability on item-level, else it is on itemset-level

    // the number of Candidate
    int NumOfCandidate = 0;

    boolean iswriteout = true;
    /**Default constructor**/
    public AlgoUUCPM() {
    }

    
    /**
     * Run the FUCPM algorithm
     * @param input the input file path
     * @param output the output file path
     * @param utilityratio minimum utility threshold ratio
     * @throws IOException exception if error while writing the file
     */
    public void runAlgorithm(String input, String output, double utilityratio) throws IOException {
        // reset MemoryLogger
        MemoryLogger.getInstance().reset();

        // input path
        this.input = input;

        // initialize the buffer for storing the current itemset
        patternBuffer = new int[BUFFERS_SIZE];

        // record the start time of the algorithm
        startTimestamp = System.currentTimeMillis();

        // create a writer object to write results to file
        writer = new BufferedWriter(new FileWriter(output));

        // for storing the current sequence number
        int NumberOfSequence = 0;

        // for storing the utility of all sequences
        double totalUtility = 0;

        BufferedReader myInput = null;
        String thisLine;

        // the database of SIL of each sequence
        List<SeqInfoList> dataset = new ArrayList<>();
        // for storing the instanceChain (composed of instanceLists) of each 1-sequence
        Map<Integer,ArrayList<InstanceList>> mapItemIC = new HashMap<>();

        //for storing the global utility and PIWU of each 1-sequence
        Map<Integer,Double> mapItemUtility = new HashMap<>();
        Map<Integer,Double> mapItemPIWU = new HashMap<>();

        //record the utility of each q-seq
        ArrayList<Double> qSeqUtility = new ArrayList();

        //record the distinct items contained by each q-seq
        HashMap<Integer, HashSet<Integer>> qSeqContainItem = new HashMap<>();
        //record the sum utility of an item in a q-seq. <key:item, value:<key:sid,value:sum utility>>
        HashMap<Integer,HashMap<Integer,Double>> mapItemSumUtility = new HashMap<>();
        //used to record the maximum probability of each item in each sequence.
        HashMap<Integer,HashMap<Integer,Double>> mapItemMaxProbability = new HashMap<>();
        //used to record the maximum probability of all items in each sequence.
        HashMap<Integer, Double> mapItemMaxSeqProbability = new HashMap<>(); 
        
        int maxlen = 0;

        /***** First scan, calculate the initial PIWU of each item, and the sum utility of each item in each q-seq *****/
        try{
            // prepare the object for reading the file
            myInput = new BufferedReader(new InputStreamReader(new FileInputStream(new File(input))));
            //q-seq id
            int Sid = 0;

            while ((thisLine = myInput.readLine()) != null) {
                // if the line is a comment, is  empty or is a kind of metadata
                if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%' || thisLine.charAt(0) == '@') {
                    continue;
                }
                // split the sequence according to the " " separator
                String tokens[] = thisLine.split(" ");
                // get the sequence utility (the last token on the line)
                String sequenceUtilityString = tokens[tokens.length-1];
                int positionColons = sequenceUtilityString.indexOf(':');
                double sequenceUtility = Double.parseDouble(sequenceUtilityString.substring(positionColons+1));
                double sequenceMaxPro = 0.0;
                //added by czf 改
                //record the utility of this q-seq
                qSeqUtility.add(sequenceUtility);
                //handle each token in this q-seq
                for(int i=0; i< tokens.length - 4; i++) {  //-2与SUtility间隔为2个空格，所以tokens.length-4
                    String currentToken = tokens[i];
                    // if empty, continue to next token
                    if (currentToken.length() == 0) {
                        continue;
                    }
                    if (currentToken.equals("-1")) {
                        continue;
                    } else {
                        // We will extract the item from the string:
                        int positionLeftBracketString1 = currentToken.indexOf('[');
                        int positionRightBracketString1 = currentToken.indexOf(']');
                        String itemString = currentToken.substring(0, positionLeftBracketString1);
                        Integer item = Integer.parseInt(itemString);
                        String utilityString = currentToken.substring(positionLeftBracketString1+1, positionRightBracketString1);
                        Double itemUtility = Double.parseDouble(utilityString);
                        int positionLeftBracketString2 = currentToken.indexOf('(');
                        int positionRightBracketString2 = currentToken.indexOf(')');
                        String ProbabilityString = currentToken.substring(positionLeftBracketString2+1, positionRightBracketString2);
                        double itemProbability = Double.parseDouble(ProbabilityString);
                        sequenceMaxPro = Math.max(sequenceMaxPro, itemProbability);
                                                
                        //add the item to qSeqContainItem
                        if(qSeqContainItem.get(Sid) == null){
                            qSeqContainItem.put(Sid, new HashSet<>());
                        }
                        qSeqContainItem.get(Sid).add(item);

                        //record the sum utility of the item in this q-seq
                        if(mapItemSumUtility.get(item) ==null){ //if it is the first occurrence of item in the database
                        	//mapItemSumUtility: the total utility of the item stored in the whole database
                                HashMap<Integer,Double> innermap = new HashMap<>();
                                innermap.put(Sid, itemUtility);
                                HashMap<Integer,Double> innermap_p = new HashMap<>();
                                innermap_p.put(Sid,itemProbability);
                                mapItemSumUtility.put(item,innermap);
                                mapItemMaxProbability.put(item, innermap_p);
//                                System.out.println(itemProbability);
                        }else{
                            HashMap<Integer, Double> innermap = mapItemSumUtility.get(item); 
                            HashMap<Integer, Double> innermap_p = mapItemMaxProbability.get(item); 
                            if(innermap.get(Sid) == null){ //if it is the first occurrence of item in this q-seq
                                innermap.put(Sid, itemUtility);
                                innermap_p.put(Sid,itemProbability);
                                //innermap_p.put(Sid, itemProbability);
//                                System.out.println(itemProbability);
                            }else{
                            	innermap.put(Sid,innermap.get(Sid) + itemUtility);
//                            	System.out.print(innermap_p.get(Sid));
                            	innermap_p.put(Sid, Math.max(innermap_p.get(Sid), itemProbability));
//                            	System.out.println(" --> " + innermap_p.get(Sid));
                            }
                        }
                    }
                }

                //update global PIWU of each item
                if(GUIRPstrategy) {
	                for(Integer item : qSeqContainItem.get(Sid)){
	                    if(mapItemPIWU.get(item)==null){
	                        mapItemPIWU.put(item, sequenceUtility * sequenceMaxPro); 
	                    }else{
	                        mapItemPIWU.put(item, sequenceUtility * sequenceMaxPro + mapItemPIWU.get(item));
	                    }
	                }
                }else {
                	for(Integer item : qSeqContainItem.get(Sid)){
	                    if(mapItemPIWU.get(item)==null){
	                        mapItemPIWU.put(item, sequenceUtility); 
	                    }else{
	                        mapItemPIWU.put(item, sequenceUtility + mapItemPIWU.get(item));
	                    }
	                }
                }
                //added by czf
                mapItemMaxSeqProbability.put(Sid, sequenceMaxPro);
                Sid++;
                //update total utility
                totalUtility+=sequenceUtility;
            }
        }catch (Exception e) {
            // catches exception if error while reading the input file
            e.printStackTrace();
        }finally {
            if(myInput != null){
                // close the input file
                myInput.close();
            }
            if(DEBUG==1) {
                for (Entry<Integer, Double> entry : mapItemPIWU.entrySet()) {
                    System.out.println("PIWU:");
                    System.out.println(entry.getKey() + " " + entry.getValue());
                }
            }
        }//First scan finished.
        //set minimum utility threshold
        minUtility = utilityratio * totalUtility;
        System.out.print("Mode: ");
        if(item_Level == true) {
        	System.out.println("Item-Level");
        }
        else {
        	System.out.println("Itemset-Level");
        }
        System.out.println("GUIRPstrategy: " + GUIRPstrategy);
        System.out.println("UEIPstrategy: " + UEIPstrategy);
        System.out.println("totalutility: "+ totalUtility);
        System.out.println("utilityratio: "+ utilityratio+" Threshold: "+minUtility);
        //the set of unpromising items
        HashSet<Integer> unpromisingItems = new HashSet<>();

        /***** Calculate Recurrent PIWU *****/
        //Recurrent PIWU: calculate the PIWU of each item recurrently.
        //The utility of unpromising items is set to zero, so the utility of each q-seq may reduce. Thus, PIWU of each item may also reduce.
        while(true){
            boolean flag = false; //record whether new unpromising item produces. true: a new unpromising item generates.
            for(Entry<Integer, Double> entry : mapItemPIWU.entrySet()){
                int item = entry.getKey();
                double itemSwu = entry.getValue();
                if(itemSwu < minUtility && unpromisingItems.contains(item)==false){
                    flag = true;
                    unpromisingItems.add(item);
                    //update utility of the q-seq which contains this unpromising item
                    for(Entry<Integer, Double> element : mapItemSumUtility.get(item).entrySet()){
                        //id of the q-seq which contains this unpromising item
                        int sid = element.getKey();
                        //sum utility of the unpromising item in sid
                        double sumUtility = element.getValue();
                        //update utility of this q-seq
                        qSeqUtility.set(sid, qSeqUtility.get(sid) - sumUtility);
                        //update PIWU of the items contained in this q-seq
                        for(int distinctItem: qSeqContainItem.get(sid)){
                            mapItemPIWU.put(distinctItem, mapItemPIWU.get(distinctItem) - sumUtility * mapItemMaxSeqProbability.get(sid));
                        }
                    }
                }
            }
            if(flag == false) break;
        }

        System.out.println("Num of unpromising items:" + unpromisingItems.size() + " Num of all distinct items:" + mapItemPIWU.size());

        if(DEBUG==2) {
            System.out.println("unpromising items:");
            for(int item : unpromisingItems){
                System.out.println(item+" ");
            }
            System.out.println("Recurrent PIWU:");
            for (Entry<Integer, Double> entry : mapItemPIWU.entrySet()) {
                System.out.println(entry.getKey() + ": " + entry.getValue());
            }
            System.out.println("q-seq utility:");
            for (int i=0; i<qSeqUtility.size(); i++) {
                System.out.println(i + ": " + qSeqUtility.get(i));
            }
        }

        /***** Second Scan, Construct the SeqInfoList, and InstanceList of promising 1-sequences (i.e. promising items) ******/
        try {
            // prepare the object for reading the file
        	//准备读取文件的对象
            myInput = new BufferedReader(new InputStreamReader(new FileInputStream(new File(input))));

            // We will read each sequence in buffers.
            // The first buffer will store the items of a sequence and the -1 between them
            int[] itemBuffer = new int[BUFFERS_SIZE];
            // The second buffer will store the utility of items in a sequence and the -1 between them
            double[] utilityBuffer = new double[BUFFERS_SIZE];
            double[] probabilityBuffer = new double[BUFFERS_SIZE];
            // The following variable will contain the length of the data stored in the two previous buffer
            int itemBufferLength;

            // for each line (q-seq) until the end of file
            while ((thisLine = myInput.readLine()) != null) {
                // if the line is  a comment, is  empty or is a kind of metadata
                if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%' || thisLine.charAt(0) == '@') {
                    continue;
                }

                // for storing the instanceList of each 1-sequence in the current q-seq
                Map<Integer, InstanceList> mapItemIL = new HashMap<>();

                //for storing utility of each 1-sequence in the current q-seq
                Map<Integer,Double> mapItemU = new HashMap<>();

                // We reset the following buffer length to zero because we are reading a new sequence.
                itemBufferLength = 0;

                // split the sequence according to the " " separator
                String tokens[] = thisLine.split(" ");

                // get the sequence utility
                double sequenceUtility = qSeqUtility.get(NumberOfSequence);

                // the number of itemsets
                int nbItemsets = 1;

                // For each token on the line except the last three tokens (the -1 -2 and SUtility).
                for(int i=0; i< tokens.length - 4; i++) {
                    String currentToken = tokens[i];
                    // if empty, continue to next token
                    if(currentToken.length() == 0) {
                        continue;
                    }
                    // if the current token is -1 ,the ending sign of an itemset
                    if(currentToken.equals("-1")) {
                        // We store the -1 in the respective buffers
                        itemBuffer[itemBufferLength] = -1;
                        utilityBuffer[itemBufferLength] = -1;
                        probabilityBuffer[itemBufferLength] = -1;
                        // We increase the length of the data stored in the buffers
                        itemBufferLength++;

                        // we update the number of itemsets in that sequence that are not empty
                        nbItemsets++;
                    }else {
                        //We need to finish the following three tasks if the current token is an item

                        /* Task 1: record the utility for constructing the SIL later */
                        // extract the item from the string:
                        int positionLeftBracketString1 = currentToken.indexOf('[');
                        int positionRightBracketString1 = currentToken.indexOf(']');
                        String itemString = currentToken.substring(0, positionLeftBracketString1);
                        Integer item = Integer.parseInt(itemString);
                        // extract the utility from the string:
                        String utilityString = currentToken.substring(positionLeftBracketString1 + 1, positionRightBracketString1);
                        Double itemUtility = Double.parseDouble(utilityString);
                        //added by czf 改
                        int positionLeftBracketString2 = currentToken.indexOf('(');
                        int positionRightBracketString2 = currentToken.indexOf(')');
                        String ProbabilityString = currentToken.substring(positionLeftBracketString2+1, positionRightBracketString2);
                        double itemProbability = Double.parseDouble(ProbabilityString);
                        //itemProbabilityMax = Math.max(itemProbabilityMax, itemProbability);
                        
                        // if the item is unpromising, we set its utility to 0
                        if (unpromisingItems.contains(item) == true) {
                            itemUtility = 0.0;
                            itemProbability = 1.0; 
                        }
                        // We store the item and its utility in the buffers for temporarily storing the sequence
                        itemBuffer[itemBufferLength] = item;
                        utilityBuffer[itemBufferLength] = itemUtility;
                        // added by czf
                        probabilityBuffer[itemBufferLength] = itemProbability;
                        itemBufferLength++;

                        // If the 1-seq is promising
                        if (itemUtility != 0) {
                            /* Task 2: Construct InstanceList of promising 1-seq */
                            // if the promising item appears in the current q-seq for the first time
                            if (mapItemIL.get(item) == null) {
                                InstanceList tempUL = new InstanceList(); 
                                tempUL.set_sid(NumberOfSequence);
                                tempUL.add(nbItemsets - 1, itemUtility, itemProbability); 
                                mapItemIL.put(item, tempUL);
                            }else {
                                InstanceList tempUL = mapItemIL.get(item);
                                tempUL.add(nbItemsets - 1, itemUtility, itemProbability); 
                                mapItemIL.put(item, tempUL);
                            } //mapItemIL: for storing the instanceList of each 1-sequence in the current q-seq

                            /* Task 3: Calculate utility of promising 1-seq */
                            // if the promising item appears in the current q-seq for the first time
                            if (mapItemU.get(item) == null) {
                            	//mapItemU.put(item, itemUtility);
                                mapItemU.put(item, itemUtility * itemProbability);
                            } else {
                            	if (itemUtility * itemProbability > mapItemU.get(item)) {
                                	//mapItemU.put(item, itemUtility);
                                    mapItemU.put(item, itemUtility * itemProbability); 
                                }
                            }
                        }
                    }
                }

                //Update global variables mapItemUtility and mapItemIC according to mapItemU and mapItemIL
                
                //mapItemU: Utility of each 1- sequence in the current q-seq;
                //mapItemIL: Used to store each 1 - sequences' instanceList in current q-seq
                //mapItemUtility: Total utility of each sequence;
                //mapItemIC: Used to store each 1 - sequences' instanceChain in current q-seq (consist of instanceList); 
                
                //update mapItemUtility, czf, 不用修改
                for(Entry<Integer,Double> entry: mapItemU.entrySet()){
                    int item = entry.getKey();
                    if(mapItemUtility.get(item) == null){
                        mapItemUtility.put(item,entry.getValue());
                    }else{
                        mapItemUtility.put(item,entry.getValue() + mapItemUtility.get(item));
                    }
                }
                
                //update mapItemIC
                for(Entry<Integer, InstanceList> entry: mapItemIL.entrySet()){
                    int item = entry.getKey();
                    ArrayList<InstanceList> tempChain = new ArrayList<InstanceList>();
                    if(mapItemIC.get(item) != null)
                        tempChain = mapItemIC.get(item);
                    tempChain.add(entry.getValue());
                    //System.out.println(entry.getValue().get_sid() + " " + entry.getValue().insList);
                    mapItemIC.put(item,tempChain);
                }

                // create the SIL for current sequence
                SeqInfoList seqinfoList = new SeqInfoList(nbItemsets); 

                // This variable will represent the position in the q-sequence
                int posBuffer = 0;
                // for each itemset 每个itemset
                for(int itemset = 0; itemset < nbItemsets; itemset++) {
                    while(posBuffer < itemBufferLength) {
                        // Get the item at the current position in the sequence
                        int item = itemBuffer[posBuffer];

                        // if it is an itemset separator, we move to next position in the sequence
                        if (item == -1) {
                            posBuffer++;
                            break;
                        }
                        // else if it is an item
                        else {
                            // get the utility of the item
                            double utility = utilityBuffer[posBuffer]; //U
                            double probability = probabilityBuffer[posBuffer];
                            // We update the rest utility by subtracting the utility of the current item
                            sequenceUtility -= utility;
                            // add the item to SIL of current q-seq if the item is promising
                            if (utility != 0) {
                                seqinfoList.registerItem(itemset, item, utility, sequenceUtility, probability);
                            }
                            posBuffer++;
                        }
                    }
                } // SIL of the current q-seq has been built. 

                // We add the SIL to the sequence database.
                dataset.add(seqinfoList);
                // if in debug mode, we print the SIL that we have just built
                if(DEBUG==3) {
                    System.out.println(seqinfoList.toString());
                    System.out.println();
                }

                // we update the number of sequences
                NumberOfSequence++;
            } // finish scaning a q-sequence each time through the loop

            // if in debug mode, we print the InstanceChain of each 1-seq
            if(DEBUG==4) {
                for (Entry<Integer,ArrayList<InstanceList>> entry: mapItemIC.entrySet()){
                    System.out.println("item:"+entry.getKey());
                    for(int i=0;i<entry.getValue().size();i++){
                        System.out.println(i+"-th InstanceList:");
                        for(int j=0;j<entry.getValue().get(i).insList.size();j++){
                            System.out.print(j+"-th element: ");
                            System.out.print("sid:"+entry.getValue().get(i).get_sid());
                            System.out.print("  tid:"+entry.getValue().get(i).insList.get(j).tid);
                            System.out.print("  acu:"+entry.getValue().get(i).insList.get(j).acu);
                            System.out.print("  pro:"+entry.getValue().get(i).insList.get(j).pro);
                        }
                        System.out.println("  End of an InstanceList");
                    }
                    System.out.println("******");
                }
            }

            // if in debug mode, we print the utility of each 1-sequence
            if(DEBUG==5) {
                System.out.println("******");
                System.out.println("Utility:");
                for (Entry<Integer,Double> entry: mapItemUtility.entrySet()){
                    System.out.println(entry.getKey()+" : "+entry.getValue());
                }
            }
        } catch (Exception e) {
            // catches exception if error while reading the input file
            e.printStackTrace();
        }finally {
            if(myInput != null){
                // close the input file
                myInput.close();
            }
        }//Second scan finished.

        // check the memory usage
        MemoryLogger.getInstance().checkMemory();

        // Mine the database recursively using the FUCPM procedure
        for(Entry<Integer, Double> entry : mapItemUtility.entrySet()){
            int item = entry.getKey();

            //GUIP pruning strategy, If PIWU is less than the minimum utility threshold, prune.
            if(mapItemPIWU.get(item) < minUtility) continue;

            patternBuffer[0]= item;
            patternBuffer[1] = -1;
            patternBuffer[2] = -2;

            NumOfCandidate++;
            //check whether the 1-seq is high-utility
            if(entry.getValue() >= minUtility){
                if(iswriteout) writeOut(patternBuffer, 1, mapItemUtility.get(item));
                patternCount++;
            }
            // recursively mine the database
           
            UUCPM(patternBuffer, 1, dataset, mapItemIC.get(item),1); 
            //PatternBuffer as prefix, 1 as prefix length, data set, mapItemIC.get(item) as prefix instanceChain, 1 as itemCount.
        }

        // check the memory usage again and close the file.
        MemoryLogger.getInstance().checkMemory();

        writer.write("=============  UUCPM ALGORITHM v1.0 - STATS ==========\n");
        writer.write(" Minimum utility threshold ~ " + minUtility + " \n");
        writer.write(" Total time ~ " + (System.currentTimeMillis() - startTimestamp)/1000 + " s\n");
        writer.write(" Max memory ~ " + MemoryLogger.getInstance().getMaxMemory() + " MB\n");
        writer.write(" Number of candidates ~ " + NumOfCandidate + " \n");
        writer.write(" Number of HUCSPs ~ " + patternCount + " \n");
        writer.write(" Number of unpromising items ~ " + unpromisingItems.size() + " \n");
        writer.write(" Number of distinct items ~ " + mapItemPIWU.size() + " \n");

        // close output file
        writer.close();
    }

    //	This inner class is used to store the information of candidates after concatenating the item.
    public class ItemConcatenation {
        // utility
        double utility;
        // projected database InstanceChain
        ArrayList<InstanceList> IChain;
        // Candidate sequence after concatenating the item
        public int[] candidate;
        // length of candidate after concatenating the item
        int candidateLength;
        // added by czf, probability
        double probability;

        // Constructor
        public ItemConcatenation(double utility, ArrayList<InstanceList> UChain, int[] candidate, int candidateLength){
            this.utility = utility;
            this.IChain = UChain;
            this.candidateLength = candidateLength;
            this.candidate = new int[BUFFERS_SIZE];
            System.arraycopy(candidate, 0, this.candidate, 0, candidateLength);
        }
 
    }

    /**
     * construct InstanceChain of candidates
     * @param Candidate candidate sequence
     * @param CandidateLength length of Candidate
     * @param database  SIL of all q-seq
     * @param instanceChain InstanceChain of the prefix of Candidate
     * @param kind 0:i-Concatenate，1:s-Concatenate
     */
    private ItemConcatenation ConstructInstanceChain(int[] Candidate, int CandidateLength, List<SeqInfoList> database, ArrayList<InstanceList> instanceChain, int kind){
        //store InstanceChain of Candidate
        ArrayList<InstanceList> ic = new ArrayList<>();
        // item: the last item of Candidate, i.e. the extension item.
        int item = Candidate[CandidateLength - 1];
//        for(int i = 0; i < CandidateLength; ++i) {
//        	System.out.print(Candidate[i] + " ");
//        }
//        System.out.print("\n");
        //the global utility of Candidate
        double Utility = 0;

        // for each InstanceList of Candidate's Prefix. Each InstanceList corresponds to one q-seq.
        for (InstanceList instanceList:instanceChain){
            // record the probabilistic utility of Candidate in current q-seq
            double LocalUtility = 0.0;
            double LocalProbability = 1.0;
            // the InstanceList of Candidate in current q-seq
            InstanceList il = new InstanceList();

            //store sid of current q-seq
            int sid = instanceList.get_sid();
            il.set_sid(sid);
            //get SIL of the current q-seq
            SeqInfoList seqinfoList = database.get(sid);
//            System.out.println("sid = " + sid);
//            System.out.println(seqinfoList.toString());

            //i-extension
            if(kind == 0){
                //construct InstanceList of Candidate in current q-seq
                //for each element in instanceList
                for (int j = 0; j < instanceList.insList.size(); j++){
                    int itemsetID = instanceList.insList.get(j).tid;
                    int itemIndex = seqinfoList.ItemsetContainItem(itemsetID, item);
                    // if the extension item appears in this itemset, we can do i-concatenate to form Candidate and get a new instance.
                    if (itemIndex != -1){
                        double PrefixUtility = instanceList.insList.get(j).acu;
//                        System.out.println("PrefixUtility: " + PrefixUtility);
                        //added and revised by czf
                        double PrefixProbability = instanceList.insList.get(j).pro;
                        il.add(itemsetID, PrefixUtility + seqinfoList.seqInfo[itemsetID].get(itemIndex).getUtility(), 
                        		PrefixProbability * seqinfoList.seqInfo[itemsetID].get(itemIndex).getProbability());
//                        System.out.println("\tadd: " + "sid: " + sid + " tid: " + itemsetID 
//                        		+ " item:" + seqinfoList.seqInfo[itemsetID].get(itemIndex).item
//                        		+ " utility:" + (PrefixUtility + seqinfoList.seqInfo[itemsetID].get(itemIndex).getUtility())
//                        		+ " " + " probability: " + PrefixProbability * seqinfoList.seqInfo[itemsetID].get(itemIndex).getProbability());
                    }
                }
                // if the current q-seq does not contain Candidate, we continue to handle the next instanceList.
                if (il.insList.size() == 0)
                    continue;
            }else{ //s-extension
                //number of itemsets in current q-seq
                int numOfItemset = seqinfoList.seqInfo.length;
                //construct InstanceList of Candidate in the current q-seq
                //for each element in instanceList
                for (int j = 0; j < instanceList.insList.size(); j++){
                    int itemsetID = instanceList.insList.get(j).tid;
                    //if it is the last itemset of current q-seq, we finish the s-extension.
                    if( itemsetID == numOfItemset - 1) break;
                    int itemIndex = seqinfoList.ItemsetContainItem(itemsetID + 1,item);
                    // if the extension item appears in the next itemset, we can do s-concatenate to form Candidate and get a new instance.
                    if (itemIndex != -1){
                        double PrefixUtility = instanceList.insList.get(j).acu;
                        //added and revised by czf
                        double PrefixProbability = instanceList.insList.get(j).pro;
                        il.add(itemsetID + 1,PrefixUtility + seqinfoList.seqInfo[itemsetID + 1].get(itemIndex).getUtility()
                        		, PrefixProbability * seqinfoList.seqInfo[itemsetID + 1].get(itemIndex).getProbability());
                    }
                }
                // if the current q-seq does not contain Candidate, we continue to handle the next instanceList.
                if (il.insList.size() == 0)
                    continue;
            }
            // calculate utility of Candidate in the current q-seq
            for (int i = 0; i < il.insList.size(); i++){
                InstanceList.InstanceElement ie = il.insList.get(i);
                if (ie.acu * ie.pro > LocalUtility)
                    LocalUtility = ie.acu * ie.pro;
            }

            //update the global Utility of Candidate
            Utility += LocalUtility;

            //add InstanceList to InstanceChain
            ic.add(il);
        }

        // if in debug mode, we print the InstanceChain that we have just built
        if(DEBUG==6){
            System.out.println("**********************");
            System.out.print("Candidate: ");
            for (int i=0;i<CandidateLength;i++){
                System.out.print(Candidate[i]+" ");
            }
            System.out.println();
            System.out.println(" global Utility:"+Utility);
            for (int i=0;i<ic.size();i++){
                System.out.println(i+"-th InstanceList:");
                for(int j=0;j<ic.get(i).insList.size();j++){
                    System.out.print("Element"+j+": ");
                    System.out.print("sid:"+ic.get(i).get_sid());
                    System.out.print("  tid:"+ic.get(i).insList.get(j).tid);
                    System.out.print("  acu:"+ic.get(i).insList.get(j).acu);
                    System.out.print("  pro:"+ic.get(i).insList.get(j).pro);
                }
                System.out.println("  End of a InstanceList");
            }
            System.out.println("#####################");
        }

        Candidate[CandidateLength] = -1;
        Candidate[CandidateLength+1] = -2;

        //return the ItemConcatenation of Candidate
        return new ItemConcatenation(Utility, ic, Candidate, CandidateLength);
    }

    /**
     * recursive pattern growth function
     * @param prefix prefix sequence
     * @param prefixLength length of prefix
     * @param database  SIL of all q-seq
     * @param instanceChain InstanceChain of prefix
     * @param itemCount number of items in prefix
     */
    private void UUCPM(int[] prefix, int prefixLength, List<SeqInfoList> database, ArrayList<InstanceList> instanceChain,
    		int itemCount) throws IOException {

        if(DEBUG==7){
                // Print the current prefix
                for(int i = 0; i < prefixLength; i++){
                    System.out.print(prefix[i] + " ");
                }
                System.out.println("TmpMinUtility:" + minUtility);
        }

        //for storing global PIEU of i-extension items and s-extension items.
        //they are also ilist and slist.
        Map<Integer,double[]> mapiItemPIEU = new HashMap<>();
        Map<Integer,double[]> mapsItemPIEU = new TreeMap<>();

        /***** Construct ilist and slist *****/ //construct ilist and slist
        //scan prefix-projected DB once to find items to be concatenated
        
        for (InstanceList instanceList : instanceChain) {
            SeqInfoList seqinfoList = database.get(instanceList.get_sid());

            //record the last item of prefix
            int item = prefix[prefixLength-1];

            // store the local PIEU of the i-extension items in current q-seq
            Map<Integer,double[]> mapiItemLocalPIEU = new HashMap<Integer, double[]>();
            // store the local PIEU of the s-extension items in current q-seq
            Map<Integer,double[]> mapsItemLocalPIEU = new HashMap<Integer, double[]>();

            /***** Construct ilist *****/
            // put i-extension items into ilist and update the global variable mapiItemPIEU
            for (int j = 0; j < instanceList.insList.size(); j++){
                int itemsetID = instanceList.insList.get(j).tid;
                //find i-extension items in current itemset
                for(int i = 0;i < seqinfoList.seqInfo[itemsetID].size(); i++){
                    //only the items whose lexicographical order is larger than that of item can be added to ilist
                    if(seqinfoList.seqInfo[itemsetID].get(i).getItem() <= item) continue;
                    int ConItem = seqinfoList.seqInfo[itemsetID].get(i).getItem();
                        //calculate PIEU of ConItem
                        double prefixUtility = instanceList.insList.get(j).acu;
                        double Probability_i = seqinfoList.seqInfo[itemsetID].get(i).getProbability();
                        double Probability_max = 0.0;
                        double currentPIEU_withoutpr = 0.0;
                        double currentPIEU = 0.0;
                        if(UEIPstrategy && item_Level) {
	                        for(int index = i + 1; index < seqinfoList.seqInfo[itemsetID].size(); index++) {
	                        	Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID].get(index).getProbability());
	                        }
	                        if(itemsetID != seqinfoList.seqInfo.length - 1) {
		                        for(int index = 0; index < seqinfoList.seqInfo[itemsetID + 1].size(); index++) {
		                        	Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID + 1].get(index).getProbability());
		                        }
	                        }
	                        currentPIEU_withoutpr = (prefixUtility + seqinfoList.seqInfo[itemsetID].get(i).getUtility() + 
	                        		seqinfoList.seqInfo[itemsetID].get(i).getRestutility())
	                        		 * instanceList.insList.get(j).pro * Probability_i;
	                        currentPIEU = currentPIEU_withoutpr * Probability_max;
                        }
                        else if(UEIPstrategy && !item_Level) {
	                        Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID].get(0).getProbability());
		                    Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID + 1].get(0).getProbability());
		                    currentPIEU_withoutpr = (prefixUtility + seqinfoList.seqInfo[itemsetID].get(i).getUtility() + 
	                        		seqinfoList.seqInfo[itemsetID].get(i).getRestutility())
	                        		 * instanceList.insList.get(j).pro * Probability_i;
		                    currentPIEU = currentPIEU_withoutpr * Probability_max;
                        }
                        else {
                        	currentPIEU_withoutpr = (prefixUtility + seqinfoList.seqInfo[itemsetID].get(i).getUtility() + 
                            		seqinfoList.seqInfo[itemsetID].get(i).getRestutility());
                        	currentPIEU = currentPIEU_withoutpr;
                        }
                        double[] currentPIEUs = new double[2];
                        currentPIEUs[0] = currentPIEU_withoutpr;
                        currentPIEUs[1] = currentPIEU;
                        //if ConItem appears in current q-seq for the first time
                        if(mapiItemLocalPIEU.get(ConItem) == null) {
                            mapiItemLocalPIEU.put(ConItem, currentPIEUs);
                            if (mapiItemPIEU.get(ConItem) == null) { //if ConItem appears in database for the first time 如果ConItem第一次出现在数据库中
                                mapiItemPIEU.put(ConItem, currentPIEUs);
                            } else {
                                double tmpPIEU_withoutpr = mapiItemPIEU.get(ConItem)[0];
                                double tmpPIEU = mapiItemPIEU.get(ConItem)[1];
                                double[] tmpPIEUs = new double[2];
                                tmpPIEUs[0] = currentPIEU_withoutpr + tmpPIEU_withoutpr;
                                tmpPIEUs[1] = currentPIEU + tmpPIEU;
                                mapiItemPIEU.put(ConItem, tmpPIEUs);
                            }
                        }else{ //if ConItem has already appeared in current q-seq.
                            if(currentPIEU > mapiItemLocalPIEU.get(ConItem)[1]){
                            	double tmpGlobalPIEU_withoutpr = mapiItemPIEU.get(ConItem)[0];
                                double tmpGlobalPIEU = mapiItemPIEU.get(ConItem)[1];
                                double[] tmpGlobalPIEUs = new double[2];
                                tmpGlobalPIEUs[0] = tmpGlobalPIEU_withoutpr - mapiItemLocalPIEU.get(ConItem)[0] + currentPIEU_withoutpr;
                                tmpGlobalPIEUs[1] = tmpGlobalPIEU - mapiItemLocalPIEU.get(ConItem)[1] + currentPIEU;
                                mapiItemPIEU.put(ConItem, tmpGlobalPIEUs);
                                mapiItemLocalPIEU.put(ConItem, currentPIEUs);
                            }
                        }

                }
            }

            /***** Construct slist *****/
            //put s-extension items into slist and get the items to be s-concatenated
            // for each element in instanceList
            for (int j = 0; j < instanceList.insList.size(); j++) {
                int itemsetID = instanceList.insList.get(j).tid;
                //if it is the last itemset of current q-seq, we finish finding s-extension items.
                if(itemsetID == seqinfoList.seqInfo.length - 1) break;
                //find s-extension items in the next itemset
                for(int i = 0; i < seqinfoList.seqInfo[itemsetID + 1].size(); i++){
                    int ConItem = seqinfoList.seqInfo[itemsetID + 1].get(i).getItem();
                        //calculate PIEU of ConItem
                        //int prefixUtility = (int)Math.round(instanceList.insList.get(j).acu * instanceList.insList.get(j).pro);
                    	double prefixUtility = instanceList.insList.get(j).acu;
                        double Probability_i = seqinfoList.seqInfo[itemsetID + 1].get(i).getProbability();
                        double Probability_max = 0.0;
                        double currentPIEU_withoutpr = 0.0;
                        double currentPIEU = 0.0;
                        
                        if(UEIPstrategy && item_Level) {
                        	if(itemsetID != seqinfoList.seqInfo.length - 2) {
		                        for(int index = i + 1; index < seqinfoList.seqInfo[itemsetID + 1].size(); index++) {
		                        	Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID + 1].get(index).getProbability());
		                        }
		                        for(int index = 0; index < seqinfoList.seqInfo[itemsetID + 2].size(); index++) {
		                        	Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID + 2].get(index).getProbability());
		                        }
                        	}
                        	currentPIEU_withoutpr = (prefixUtility + seqinfoList.seqInfo[itemsetID + 1].get(i).getUtility() + 
                            		seqinfoList.seqInfo[itemsetID + 1].get(i).getRestutility())
                            		* instanceList.insList.get(j).pro * Probability_i;
                            currentPIEU = currentPIEU_withoutpr  * Probability_max;
                        }
                        else if(UEIPstrategy && !item_Level) {
                        	if(itemsetID != seqinfoList.seqInfo.length - 2) {
		                        Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID + 1].get(0).getProbability());
		                        Probability_max = Math.max(Probability_max, seqinfoList.seqInfo[itemsetID + 2].get(0).getProbability());
                        	}
                        	currentPIEU_withoutpr = (prefixUtility + seqinfoList.seqInfo[itemsetID + 1].get(i).getUtility() + 
                            		seqinfoList.seqInfo[itemsetID + 1].get(i).getRestutility())
                            		* instanceList.insList.get(j).pro * Probability_i;
                            currentPIEU = currentPIEU_withoutpr  * Probability_max;
                        }
                        else {
                        	currentPIEU_withoutpr = (prefixUtility + seqinfoList.seqInfo[itemsetID + 1].get(i).getUtility() + 
                            		seqinfoList.seqInfo[itemsetID + 1].get(i).getRestutility());
                            currentPIEU = currentPIEU_withoutpr;
                        }
                        double[] currentPIEUs = new double[2];
                        currentPIEUs[0] = currentPIEU_withoutpr;
                        currentPIEUs[1] = currentPIEU;
                        //if ConItem appears in current q-seq for the first time
                        if(mapsItemLocalPIEU.get(ConItem) == null) {
                            mapsItemLocalPIEU.put(ConItem, currentPIEUs);
                            if (mapsItemPIEU.get(ConItem) == null) { //if ConItem appears in database for the first time
                                mapsItemPIEU.put(ConItem, currentPIEUs);
                            } else {
                            	double tmpPIEU_withoutpr = mapsItemPIEU.get(ConItem)[0];
                                double tmpPIEU = mapsItemPIEU.get(ConItem)[1];
                                double[] tmpPIEUs = new double[2];
                                tmpPIEUs[0] = currentPIEU_withoutpr + tmpPIEU_withoutpr;
                                tmpPIEUs[1] = currentPIEU + tmpPIEU;
                                mapsItemPIEU.put(ConItem, tmpPIEUs);
                            }
                            //mapsItemLocalPIEU.put(ConItem, 1);
                        }else{ //if ConItem has already appeared in current q-seq.
                            //choose the greater PIEU
                        	if(currentPIEU > mapsItemLocalPIEU.get(ConItem)[1]){
                            	double tmpGlobalPIEU_withoutpr = mapsItemPIEU.get(ConItem)[0];
                                double tmpGlobalPIEU = mapsItemPIEU.get(ConItem)[1];
                                double[] tmpGlobalPIEUs = new double[2];
                                tmpGlobalPIEUs[0] = tmpGlobalPIEU_withoutpr - mapsItemLocalPIEU.get(ConItem)[0] + currentPIEU_withoutpr;
                                tmpGlobalPIEUs[1] = tmpGlobalPIEU - mapsItemLocalPIEU.get(ConItem)[1] + currentPIEU;
                                mapsItemPIEU.put(ConItem, tmpGlobalPIEUs);
                                mapsItemLocalPIEU.put(ConItem, currentPIEUs);
                            }
                        }
                }
            }
        }//Finish constructing ilist and slist.

        // for temporarily storing information of candidates after extension
        ItemConcatenation ItemCom;

        /***** I-Extension *****/ 
        // perform I-Extension to grow the pattern larger.
        for (Entry<Integer,double[]> entry : mapiItemPIEU.entrySet()){
            int item = entry.getKey();
            double ieu_withoutpr = entry.getValue()[0];
            double ieu = entry.getValue()[1];

            // PIEU pruning strategy
            if (ieu_withoutpr < minUtility){
                continue;
            }//Actually, it is only valid when the candidate pattern length is 1.

            //construct the candidate after extension
            prefix[prefixLength] = item;

            //construct InstanceChain of the candidate
            if (itemCount + 1 <= maxPatternLength){
                ItemCom = ConstructInstanceChain(prefix, prefixLength + 1, database, instanceChain, 0);
                NumOfCandidate++;
//                System.out.print("candidate: [");
//                for (int i = 0; i < ItemCom.candidateLength; i++) {
//                	System.out.print((prefix[i] + " "));
//                }
//                System.out.println("-1]");
                //check whether the candidate is high-utility
                if(ItemCom.utility >= minUtility){                	
                    if(iswriteout) writeOut(ItemCom.candidate, ItemCom.candidateLength, ItemCom.utility);
                    patternCount++;
                }
                // PIEU pruning strategy
                if (ieu < minUtility){
                    continue;
                }

                //mine the database recursively using the FUCPM procedure
                UUCPM(ItemCom.candidate, ItemCom.candidateLength, database, ItemCom.IChain, itemCount+1);
            }
        }

        /***** S-Extension *****/
        // perform S-Extension to grow the pattern larger.
        for (Entry<Integer, double[]> entry : mapsItemPIEU.entrySet()){
            int item = entry.getKey();
            double ieu_withoutpr = entry.getValue()[0];
            double ieu = entry.getValue()[1];

            //PIEU pruning strategy
            if (ieu_withoutpr < minUtility){
                continue;
            }//Actually, it is only valid when the candidate pattern length is 1.

            //construct the candidate after extension
            prefix[prefixLength] = -1;
            prefix[prefixLength+1] = item;

            //construct InstanceChain of the candidate
            if (itemCount + 1 <= maxPatternLength){
                ItemCom = ConstructInstanceChain(prefix,prefixLength + 2, database, instanceChain,1);
                NumOfCandidate++;
//                System.out.print("candidate: [");
//                for (int i = 0; i < ItemCom.candidateLength; i++) {
//                	System.out.print((prefix[i] + " "));
//                }
//                System.out.println("-1]");
                //check whether the candidate is high-utility
                if(ItemCom.utility >= minUtility){
                    if(iswriteout) writeOut(ItemCom.candidate, ItemCom.candidateLength, ItemCom.utility);
                    patternCount++;
                }
                
                if (ieu_withoutpr < minUtility){
                    continue;
                }
                
                //mine the database recursively using the FUCPM procedure
                UUCPM(ItemCom.candidate, ItemCom.candidateLength, database, ItemCom.IChain, itemCount+1);
            }
        }

        //check the memory usage
        MemoryLogger.getInstance().checkMemory();
    }

    /**
     * Set the maximum pattern length
     * @param maxPatternLength the maximum pattern length
     */
    public void setMaxPatternLength(int maxPatternLength) {
        this.maxPatternLength = maxPatternLength;
    }

    /**
     * Method to write a high utility itemset to the output file.
     //* @param the prefix to be written o the output file
     * @param utility the utility of the prefix concatenated with the item
     * @param prefixLength the prefix length
     */
    private void writeOut(int[] prefix, int prefixLength,  double utility) throws IOException {
        // increase the number of high utility itemsets found
        //patternCount++;

        StringBuilder buffer = new StringBuilder();

        // If the user wants to save in SPMF format
        if(SAVE_RESULT_EASIER_TO_READ_FORMAT == false) {
            // append each item of the pattern
            for (int i = 0; i < prefixLength; i++) {
                buffer.append(prefix[i]);
                buffer.append(' ');
            }

            // append the end of itemset symbol (-1) and end of sequence symbol (-2)
            buffer.append("-1 #UTIL: ");
            // append the utility of the pattern
            buffer.append(String.format("%.2f", utility));
        }
        else {
            // Otherwise, if the user wants to save in a format that is easier to read for debugging.
            // Append each item of the pattern
            buffer.append('<');
            buffer.append('(');
            for (int i = 0; i < prefixLength; i++) {
                if(prefix[i] == -1) {
                    buffer.append(")(");
                }else {
                    buffer.append(prefix[i]);
                }
            }
            buffer.append(")>:");
            buffer.append(String.format("%.2f", utility));
        }

        // output, czf
//        System.out.println(buffer.toString());

        
        // write the pattern to the output file
        writer.write(buffer.toString());
        writer.newLine();

        // if in debugging mode, then also print the pattern to the console
        if(DEBUG==8) {
            System.out.println(" SAVING : " + buffer.toString());
            System.out.println();

            // check if the calculated utility is correct by reading the file
            // for debugging purpose
            checkIfUtilityOfPatternIsCorrect(prefix, prefixLength, utility);
        }
    }

    /**
     * This method check if the utility of a pattern has been correctly calculated for
     * debugging purposes. It is not designed to be efficient since it is just used for
     * debugging.
     * @param prefix a pattern stored in a buffer
     * @param prefixLength the pattern length
     * @param utility the utility of the pattern
     * @throws IOException if error while writting to file
     */
    private void checkIfUtilityOfPatternIsCorrect(int[] prefix, int prefixLength, double utility) throws IOException {
        double calculatedUtility = 0;

        BufferedReader myInput = new BufferedReader(new InputStreamReader( new FileInputStream(new File(input))));
        // we will read the database
        try {
            // prepare the object for reading the file

            String thisLine;
            // for each line (q-seq) until the end of file
            while ((thisLine = myInput.readLine()) != null) {
                // if the line is  a comment, is  empty or is a kind of metadata
                if (thisLine.isEmpty() == true || thisLine.charAt(0) == '#' || thisLine.charAt(0) == '%' || thisLine.charAt(0) == '@') {
                    continue;
                }

                // split the sequence according to the " " separator
                String tokens[] = thisLine.split(" ");

                int tokensLength = tokens.length -3;

                int[] sequence = new int[tokensLength];
                double[] sequenceUtility = new double[tokensLength];

                // Copy the current sequence in the sequence buffer.
                // For each token on the line except the last three tokens
                // (the -1 -2 and sequence utility).
                for(int i=0; i< tokensLength; i++) {
                    String currentToken = tokens[i];

                    // if empty, continue to next token
                    if(currentToken.length() == 0) {
                        continue;
                    }

                    // read the current item
                    int item;
                    double itemUtility;

                    // if the current token is -1
                    if(currentToken.equals("-1")) {
                        item = -1;
                        itemUtility = 0;
                    }else {
                        // if  the current token is an item
                        //  We will extract the item from the string:
                        int positionLeftBracketString = currentToken.indexOf('[');
                        int positionRightBracketString = currentToken.indexOf(']');
                        String itemString = currentToken.substring(0, positionLeftBracketString);
                        item = Integer.parseInt(itemString);

                        // We also extract the utility from the string:
                        String utilityString = currentToken.substring(positionLeftBracketString+1, positionRightBracketString);
                        itemUtility = Integer.parseInt(utilityString);
                    }
                    sequence[i] = item;
                    sequenceUtility[i] = itemUtility;
                }

                // For each position of the sequence
                double util = tryToMatch(sequence,sequenceUtility, prefix, prefixLength, 0, 0, 0);
                calculatedUtility += util;
            }
        } catch (Exception e) {
            // catches exception if error while reading the input file
            e.printStackTrace();
        }finally {
            if(myInput != null){
                // close the input file
                myInput.close();
            }
        }

        if(calculatedUtility != utility) {
            System.out.print(" ERROR, WRONG UTILITY FOR PATTERN : ");
            for(int i=0; i<prefixLength; i++) {
                System.out.print(prefix[i]);
            }
            System.out.println(" utility is: " + utility + " but should be: " + calculatedUtility);
            System.in.read();
        }
    }

    /**
     * This is some code for verifying that the utility of a pattern is correctly calculated
     * for debugging only. It is not efficient. But it is a mean to verify that
     * the result is correct.
     * @param sequence a sequence (the items and -1)
     * @param sequenceUtility a sequence (the utility values and -1)
     * @param prefix the current pattern stored in a buffer
     * @param prefixLength the current pattern length
     * @param prefixPos the position in the current pattern that we will try to match with the sequence
     * @param seqPos the position in the sequence that we will try to match with the pattenr
     * @param utility the calculated utility until now
     * @return the utility of the pattern
     */
    private double tryToMatch(int[] sequence, double[] sequenceUtility, int[] prefix,	int prefixLength,
                           int prefixPos, int seqPos, double utility) {

        // Note: I do not put much comment in this method because it is just
        // used for debugging.

        List<Double> otherUtilityValues = new ArrayList<Double>();

        // try to match the current itemset of prefix
        int posP = prefixPos;
        int posS = seqPos;

        int previousPrefixPos = prefixPos;
        double itemsetUtility = 0;
        while(posP < prefixLength & posS < sequence.length) {
            if(prefix[posP] == -1 && sequence[posS] == -1) {
                posS++;

                // try to skip the itemset in prefix
                double otherUtility = tryToMatch(sequence, sequenceUtility, prefix, prefixLength, previousPrefixPos, posS, utility);
                otherUtilityValues.add(otherUtility);

                posP++;
                utility += itemsetUtility;
                itemsetUtility = 0;
                previousPrefixPos = posP;
            }else if(prefix[posP] == -1) {
                // move to next itemset of sequence
                while(posS < sequence.length && sequence[posS] != -1){
                    posS++;
                }

                // try to skip the itemset in prefix
                double otherUtility = tryToMatch(sequence, sequenceUtility, prefix, prefixLength, previousPrefixPos, posS, utility);
                otherUtilityValues.add(otherUtility);

                utility += itemsetUtility;
                itemsetUtility = 0;
                previousPrefixPos = posP;

            }else if(sequence[posS] == -1) {
                posP = previousPrefixPos;
                itemsetUtility = 0;
                posS++;
            }else if(prefix[posP] == sequence[posS]) {
                posP++;
                itemsetUtility += sequenceUtility[posS];
                posS++;
                if(posP == prefixLength) {

                    // try to skip the itemset in prefix
                    // move to next itemset of sequence
                    while(posS < sequence.length && sequence[posS] != -1){
                        posS++;
                    }
                    double otherUtility = tryToMatch(sequence, sequenceUtility, prefix, prefixLength, previousPrefixPos, posS, utility);
                    otherUtilityValues.add(otherUtility);


                    utility += itemsetUtility;
                }
            }else if(prefix[posP] != sequence[posS]) {
                posS++;
            }
        }

        double max = 0;
        if(posP == prefixLength) {
            max = utility;
        }
        for(Double utilValue : otherUtilityValues) {
            if(utilValue > utility) {
                max = utilValue;
            }
        }
        return max;
    }

    /**
     * Print statistics about the latest execution to System.out.
     */
    public void printStatistics() {
        System.out.println("=============  UUCPM ALGORITHM v3.0 - STATS  ==========");
        System.out.println(" Minimum utility threshold ~ " + String.format("%.2f", minUtility));
        System.out.println(" Total time ~ " + (System.currentTimeMillis() - startTimestamp)/1000 + " s");
        System.out.println(" Max memory ~ " + MemoryLogger.getInstance().getMaxMemory() + " MB");
        System.out.println(" Number of candidates ~ " + NumOfCandidate);
        System.out.println(" Number of UUCSPs ~ " + patternCount);
        System.out.println("========================================================");
    }
}
