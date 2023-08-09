import java.io.IOException;

/**
 * Copyright (C), 2022-2023, Jinan University
 * FileName: MainTestUUCPM
 * Author:   Zefeng Chen
 * Date:     2023/5/20 21:20
 * Test the UUCPM algorithm.
 * 
 */
public class MainTestUUCPM {

    public static void main(String [] arg) throws IOException{

        //datasets
        String[] datasets = {"testpaper"};

        //threshold
        double [][]minutil = {
                {4e-6},
        };

        int index = 0;
        // run the algorithm
        for(String s: datasets) {
            // the input database
            String input = "input/" + s + ".txt";
            for (double i : minutil[index]) {
                //AlgoFUCPM algo = new AlgoFUCPM();
                AlgoUUCPM algo = new AlgoUUCPM();
                // set the maximum pattern length (optional)
                algo.setMaxPatternLength(1000);
                // the path for saving the patterns found
                String output = "output/UUCPM"+s+"_" + i + ".txt";
                algo.runAlgorithm(input, output, i);
                // print statistics
                algo.printStatistics();
            }
            index++;
        }
    }
}
