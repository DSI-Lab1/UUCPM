/**
 * Copyright (C), 2015-2020, HITSZ
 * FileName: InstanceList
 * Author:   qj Dai
 * Date:     2021/1/24 10:49
 * Description: InstanceList is used to store the utility and position of instances of a pattern
 */


import java.util.ArrayList;

public class InstanceList {
	
	

    public class InstanceElement {

        //the id of the itemset where the last item of the instance locates
    	//实例的最后一项所在的项集的id
        public int tid;

        //utility of the instance
        //实例的效用值
        public double acu;
        
        //added by czf
        //实例的概率
        public double pro;

        public InstanceElement(int tid, double acu) {
            this.tid = tid;
            this.acu = acu;
        }
        
        //added by czf
        public InstanceElement(int tid, double acu, double pro) {
            this.tid = tid;
            this.acu = acu;
            this.pro = pro;
        }
    }

    //InstanceList:
    ArrayList<InstanceElement> insList = new ArrayList<>(2);

    //q-sequence id
    public int sid;

    public InstanceList() {
    }

    public void add(int tid, double acu) {
        this.insList.add(new InstanceElement(tid, acu));
    }
    
    //added by czf
    public void add(int tid, double acu, double pro) {
        this.insList.add(new InstanceElement(tid, acu, pro));
    }

    public void set_sid(int sid) {
        this.sid = sid;
    }

    public int get_sid() {
        return this.sid;
    }

}