/*******************************************************************************
** ProbabilityState.java
** Part of the Java Mutual Information toolbox
** 
** Author: Adam Pocock
** Created: 20/1/2012
**
**  Copyright 2012-2016 Adam Pocock, The University Of Manchester
**  www.cs.manchester.ac.uk
**
**  This file is part of JavaMI.
**
**  JavaMI is free software: you can redistribute it and/or modify
**  it under the terms of the GNU Lesser General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  JavaMI is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public License
**  along with JavaMI.  If not, see <http://www.gnu.org/licenses/>.
**
*******************************************************************************/

package JavaMI;

import java.util.HashMap;
import java.util.Map.Entry;

/**
 * Calculates the probabilities of each state in a random variable.
 * Provides the base for all functions of one variable. Additional functions
 * include the normaliseArrays function which converts all inputs so they start
 * at 0, and the mergeArrays function which creates an array of the joint state of
 * the two input arrays.
 *
 * @author apocock
 */
public class ProbabilityState
{
  public final HashMap<Integer,Double> probMap;
  public final int maxState;
  
  /**
   * Constructor for the ProbabilityState class. Takes a data vector and calculates
   * the marginal probability of each state, storing each state/probability pair in a HashMap.
   *
   * @param  dataVector  Input vector. It is discretised to the floor of each value.
   */
  public ProbabilityState(double[] dataVector)
  {
    probMap = new HashMap<Integer,Double>();
    int vectorLength = dataVector.length;
    double doubleLength = dataVector.length;

    //round input to integers
    int[] normalisedVector = new int[vectorLength];
    maxState = normaliseArray(dataVector,normalisedVector);
   
    HashMap<Integer,Integer> countMap = new HashMap<Integer,Integer>();

    for (int value : normalisedVector)
    {
        Integer tmpKey = value;
        Integer tmpValue = countMap.remove(tmpKey);
        if (tmpValue == null)
        {
            countMap.put(tmpKey,1);
        }
        else
        {
            countMap.put(tmpKey,tmpValue + 1);
        }
    }

    for (Entry<Integer,Integer> e : countMap.entrySet())
    {
        probMap.put(e.getKey(),e.getValue() / doubleLength);
    }
  }//constructor(double[])


  /** 
   * Takes an input vector and writes an output vector
   * which is a normalised version of the input, and returns the maximum state.
   * A normalised array has min value = 0, max value = old max value - min value
   * and all values are integers
   *
   * The length of the vectors must be the same, and outputVector must be 
   * instantiated before calling this function.
   * @param inputVector The vector to normalise.
   * IMPORTANTE NOTE:
   * Hanchuan Peng's MutualInfo 0.9 MATLAB library makes the following extra logic for normalization.
   * From his code:
   * "Originally I added 0.5 before rounding, however seems the negative numbers and 
   *      positive numbers are all rounded towarded 0; hence int(-1+0.5)=0 and int(1+0.5)=1;
   *      This is unwanted because I need the above to be -1 and 1.
   *      so, I just round with 0.5 adjustment for positive and negative differently"
   * For this reason, data will be adjusted around -0.5 and 0.5 if input is positive or negative, respectively.
   * 
   * @param outputVector The normalised vector. Must be instantiated to length inputVector.length.
   * @return The maximum state from the normalised vector.
   */
  public static int normaliseArray(double[] inputVector, int[] outputVector)
  {
    int minVal = 0;
    int maxVal = 0;
    int currentValue = 0;
    int vectorLength = inputVector.length;
    
    if (vectorLength > 0)
    {
    	// Hanchuan Peng extra logic.
        if (inputVector[0] > 0 ) {
    		minVal = maxVal = (int) (inputVector[0] + 0.5);
    	} else {
    		minVal = maxVal = (int) (inputVector[0] - 0.5);
    	}
    
      for (int i = 0; i < vectorLength; i++)
      {
    	// Hanchuan Peng extra logic.
    	  if (inputVector[i] > 0 ) {
    		  currentValue = (int) (inputVector[i] + 0.5);
    	  } else {
    		  currentValue = (int) (inputVector[i] - 0.5);
    	  }
        outputVector[i] = currentValue;
        
        if (currentValue < minVal)
        {
          minVal = currentValue;
        }
        
        if (currentValue > maxVal)
        {
          maxVal = currentValue;
        }
      }/*for loop over vector*/
      
      for (int i = 0; i < vectorLength; i++)
      {
        outputVector[i] = outputVector[i] - minVal;
      }

      maxVal = (maxVal - minVal) + 1;
    }

    return maxVal;
  }//normaliseArray(double[],double[])

  /**
   * Takes in two arrays and writes the joint state of those arrays
   * to the output vector, returning the maximum joint state.
   *
   * The length of all vectors must be equal to firstVector.length
   * outputVector must be instantiated before calling this function.
   * @param firstVector The first vector.
   * @param secondVector The second vector.
   * @param outputVector The merged vector. Must be instantiated to length inputVector.length.
   * @return The maximum state from the merged vector.
   */
  public static int mergeArrays(double[] firstVector, double[] secondVector, double[] outputVector)
  {
    int[] firstNormalisedVector;
    int[] secondNormalisedVector;
    int firstNumStates;
    int secondNumStates;
    int[] stateMap;
    int stateCount;
    int curIndex;
    int vectorLength = firstVector.length;
    
    firstNormalisedVector = new int[vectorLength];
    secondNormalisedVector = new int[vectorLength];

    firstNumStates = normaliseArray(firstVector,firstNormalisedVector);
    secondNumStates = normaliseArray(secondVector,secondNormalisedVector);
    
    stateMap = new int[firstNumStates*secondNumStates];
    stateCount = 1;
    for (int i = 0; i < vectorLength; i++)
    {
      curIndex = firstNormalisedVector[i] + (secondNormalisedVector[i] * firstNumStates);
      if (stateMap[curIndex] == 0)
      {
        stateMap[curIndex] = stateCount;
        stateCount++;
      }
      outputVector[i] = stateMap[curIndex];
    }
      
    return stateCount;
  }//mergeArrays(double[],double[],double[])
  
  /**
   * A helper function which prints out any given int vector.
   * Mainly used to help debug the rest of the toolbox.
   *
   * @param vector The vector to print out.
   */
  public static void printIntVector(int[] vector)
  {
    for (int i = 0; i < vector.length; i++)
    {
      if (vector[i] > 0) {
        System.out.println("Val at i=" + i + ", is " + vector[i]);
      }
    }//for number of items in vector
  }//printIntVector(doublei[])

  /**
   * A helper function which prints out any given double vector.
   * Mainly used to help debug the rest of the toolbox.
   *
   * @param vector The vector to print out.
   */
  public static void printDoubleVector(double[] vector)
  {
    for (int i = 0; i < vector.length; i++)
    {
      if (vector[i] > 0) {
        System.out.println("Val at i=" + i + ", is " + vector[i]);
      }
    }//for number of items in vector
  }//printDoubleVector(doublei[])
}//class ProbabilityState
