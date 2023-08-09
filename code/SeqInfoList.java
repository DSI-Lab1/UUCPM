/**
 * Copyright (C), 2015-2020, HITSZ
 * FileName: SeqInfoList
 * Author:   qj Dai
 * Date:     2021/1/24 10:49
 * Description: SeqInfoList (SIL) is used to store utility and rest utility of each item in q-sequences.
 */


import java.util.ArrayList;

public class SeqInfoList {

    //store utility and rest utility of an item in the q-sequence
	//存放q-seq中一个item的效用和剩余效用
    class itemInfo{
        int item;
        double utility;
        double restUtility;
        // added by czf
        double probability;

        public itemInfo(int item, double utility, double restutility){
            this.item = item;
            this.utility = utility;
            this.restUtility = restutility;
        }
        
        public itemInfo(int item, double utility, double restutility, double probability){
            this.item = item;
            this.utility = utility;
            this.restUtility = restutility;
            this.probability = probability;
        }
        public int getItem(){
            return item;
        }
        public double getUtility(){
            return utility;
        }
        public double getRestutility(){
            return restUtility;
        }
        public double getProbability(){
            return probability;
        }
    }

    // Each element of seqInfo corresponds to an itemset of the q-sequence.
    // seqInfo的每个元素都对应于q-seq的一个项集
    ArrayList<itemInfo>[] seqInfo;

    //构造方法: 将SeqInfo的所有itemset传进去
    public SeqInfoList(int NumOfItemset){
        seqInfo = new ArrayList[NumOfItemset];
        for(int i = 0; i < NumOfItemset; i++){
            seqInfo[i] = new ArrayList<itemInfo>(1);
        }
    }

    //find the position of item in itemset
    //查找项目在项目集中的位置
    public int ItemsetContainItem(int itemsetID, int item){
        for(int i = 0; i < seqInfo[itemsetID].size(); i++){
            if(seqInfo[itemsetID].get(i).item == item) return i;
            else if(seqInfo[itemsetID].get(i).item > item) break;
        }
        return -1;
    }

    //register an item
    //注册一个项目
    public void registerItem(int itemsetID, int item, double utility, double restutility) {
        seqInfo[itemsetID].add(new itemInfo(item,utility,restutility));
    }
    
    //added by czf
    public void registerItem(int itemsetID, int item, double utility, double restutility, double probability) {
        seqInfo[itemsetID].add(new itemInfo(item, utility, restutility, probability));
    }

    //transform the SIL into string
    //把SIL转换成字符串
    public String toString() {
        StringBuffer buffer = new StringBuffer();
        buffer.append(" SIL \n");
        for(int i = 0; i < seqInfo.length; i++) {
            for(int j = 0;j < seqInfo[i].size(); j++){
                buffer.append( "item: " + seqInfo[i].get(j).item + " utility: " + seqInfo[i].get(j).utility +
                        " restutility: " + seqInfo[i].get(j).restUtility + " probability: " + seqInfo[i].get(j).probability + "   ");
            }
            buffer.append("\n");
        }
        buffer.append("\n");
        return buffer.toString();
    }
}
