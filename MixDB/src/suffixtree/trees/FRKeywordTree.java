package suffixtree.trees;

import java.util.ArrayList;
import java.util.Collections;

import sequences.MassSequence;
import suffixtree.Matching;
import suffixtree.Modification;
import suffixtree.Mutation;
import suffixtree.edges.Edge;
import suffixtree.matches.ModMatchObject;
import suffixtree.matches.MutMatchObject;
import suffixtree.misc.ProgressMeter;
import suffixtree.nodes.Node;


/**
 * FRKeywordTree stands for Forward-Reverse Keyword Tree and the main reason of
 * its existence is to speed up the mutation searches by splitting the search
 * with a forward and reverse tree. In the forward tree, mutations are not
 * allowed in the first half, while in the reverse tree mutations are not allowed
 * in the second half. Once the results are joined, we should get all the mutated
 * matches. 
 * @author jung
 *
 */
public class FRKeywordTree {

  private ArrayList<ArrayList<Integer>> queries;

  private int halfLength;
  
  
  public FRKeywordTree(ArrayList<ArrayList<Integer>> queries) {
    this.queries = queries;
  }
  
  /**
   * Get the query with the given index
   * @param index the index of the query to retrieve
   * @return the query object
   */
  public ArrayList<Integer> getQueryAt(int index) {
    return this.queries.get(index);
  }
  
  
  
  /**
   * Match this tree against the text db allowing for 1 mutation
   * @param db the database used as the text
   * @param results the array of mutated matches
   */
  public void mutationMatch(MassSequence db, ArrayList<MutMatchObject> results) {
    
    // reverse the queries and construct the reverse keywordtree
    this.halfLength = Integer.MAX_VALUE;
    ArrayList<ArrayList<Integer>> rQueries = new ArrayList<ArrayList<Integer>>(); 
    for (ArrayList<Integer> query : this.queries) {
      ArrayList<Integer> a = new ArrayList<Integer>(query); 
      Collections.reverse(a);
      rQueries.add(a);
      
      if (a.size()<this.halfLength) this.halfLength = a.size();
    }
    this.halfLength = this.halfLength / 2;

    
    
    // reverse searching
    KeywordTree tree = new KeywordTree(rQueries);
    int rootMaxEdge = tree.getRoot().getMaximumEdge().getLabel();
    int rootMinEdge = tree.getRoot().getMinimumEdge().getLabel();
    
    // From the performance standpoint, it is much faster to have a database
    // buffered and read chunks of the db at a time, and proceed to shift
    // this integer array
    ArrayList<Integer> masses = new ArrayList<Integer>();       // database buffer
    ArrayList<Character> letters = new ArrayList<Character>();  // amino acid letters
    
    int edgeCount = 0;
    long position = db.getSize()-1;
    
    ProgressMeter pm = new ProgressMeter("Seaching first pass " + db.getSize(), db.getSize(), System.out);
    while(position>=0) {
    
      masses.clear(); letters.clear();
      
      // this is the coordinate of the first position of dbb in the db 
      long dbStart = position;    
      while (db.hasMass(position)) {
        masses.add(db.getIntegerMass(position)); // items here are reversed
        letters.add(db.getCharAt(position));
        position--;
      }
      
      if (masses.size()>0) {    // only process when we have something in the buffer
        Integer[] dbb = masses.toArray(new Integer[masses.size()]);
        Character[] cdbb = letters.toArray(new Character[letters.size()]);
        
        for (int dbbStart=0; dbbStart<dbb.length; dbbStart++) {
          int firstMass=0;
          for (int extend=dbbStart; extend<dbb.length; extend++) {
            firstMass += dbb[extend];
            
            // optimization stuff
            if (firstMass < rootMinEdge) continue;
            if (firstMass > rootMaxEdge) break;
            
            int matchIndex = tree.getRoot().search(firstMass);
            if (matchIndex>=0) {
              edgeCount += this.mutationMatchR(db, dbStart, dbbStart, cdbb, dbb, extend+1, tree.getRoot().getEdgeAt(matchIndex), 1, 1, results);
            }
          }
        }
        pm.update(db.getSize()-position);  
      }
     
      position--;  // skip the blank
    }
    
    System.out.println("\nAverage number of navigated edges per position " + edgeCount/db.getSize());
    tree = null;   // clear the tree
    
    
    
    // forward searching
    tree = new KeywordTree(this.queries);
    rootMaxEdge = tree.getRoot().getMaximumEdge().getLabel();
    rootMinEdge = tree.getRoot().getMinimumEdge().getLabel();
    
    edgeCount = 0;
    position = 0;
    
    pm = new ProgressMeter("Seaching second pass " + db.getSize(), db.getSize(), System.out);
    while (position<db.getSize()) {
      
      masses.clear(); letters.clear();
      
      long dbStart = position;    // register the current beginning index in db
      while (db.hasMass(position)) {
        masses.add(db.getIntegerMass(position));
        letters.add(db.getCharAt(position));
        position++;
      }
      
      if (masses.size()>0) {
        Integer[] dbb = masses.toArray(new Integer[masses.size()]);
        Character[] cdbb = letters.toArray(new Character[letters.size()]);
      
        for (int dbbStart=0; dbbStart<dbb.length; dbbStart++) { // try different starts
          int firstMass=0;   // try to match the first mass
          for (int extend=dbbStart; extend<dbb.length; extend++) {
            firstMass += dbb[extend];
            
            // optimization stuff
            if (firstMass < rootMinEdge) continue;
            if (firstMass > rootMaxEdge) break;
            
            int matchIndex = tree.getRoot().search(firstMass);
            if (matchIndex>=0) {
              edgeCount += this.mutationMatch(db, dbStart, dbbStart, cdbb, dbb, extend+1, tree.getRoot().getEdgeAt(matchIndex), 1, 1, results);
            }
          }
        }
        pm.update(position);
      }
      
      position++;  // skip the blank
    } 
    System.out.println("\nAverage number of navigated edges per position " + edgeCount/db.getSize());
    tree=null;
  }
  
  
  

  public void blindMatch(MassSequence db, ArrayList<ModMatchObject> results) {
    
    // reverse the queries and construct the reverse keywordtree
    this.halfLength = Integer.MAX_VALUE;
    ArrayList<ArrayList<Integer>> rQueries = new ArrayList<ArrayList<Integer>>(); 
    for (ArrayList<Integer> query : this.queries) {
      ArrayList<Integer> a = new ArrayList<Integer>(query); 
      Collections.reverse(a);
      rQueries.add(a);
      
      if (a.size()<this.halfLength) this.halfLength = a.size();
    }
    this.halfLength = this.halfLength / 2;
    
    // It is faster to buffer sections of the database into buffers (arrays)
    ArrayList<Integer> masses = new ArrayList<Integer>();
    
    // reverse searching
    KeywordTree tree = new KeywordTree(rQueries);
    int rootMaxEdge = tree.getRoot().getMaximumEdge().getLabel();
    int rootMinEdge = tree.getRoot().getMinimumEdge().getLabel();
    
    String msg = String.format("Searching first pass db %d with 1 blind mod: ", db.getSize());
    ProgressMeter pm = new ProgressMeter(msg, db.getSize(), System.out);
    
    long edgeCount = 0, position = db.getSize()-1;
    
    while (position >= 0) {
      masses.clear();
      
      long dbStart = position;
      while (db.hasMass(position)) {
        masses.add(db.getIntegerMass(position));
        position--;
      }
      
      if (masses.size()>0) {
        Integer[] dbb = masses.toArray(new Integer[masses.size()]);
        
        for (int dbbStart=0; dbbStart<dbb.length; dbbStart++) {
          int firstMass = 0;
          for (int extend=dbbStart; extend<dbb.length; extend++) {
            firstMass += dbb[extend];
            
            // optimization
            if (firstMass < rootMinEdge) continue;
            if (firstMass > rootMaxEdge) break;
            
            int matchIndex = tree.getRoot().search(firstMass);
            if (matchIndex>=0) {
              edgeCount += this.modificationMatchR(db, dbStart, dbbStart, dbb, extend+1, tree.getRoot().getEdgeAt(matchIndex), 1, 1, results);             
            }
          }
          
        }
        pm.update(db.getSize()-position);
      }
      
      position--;
    }
    System.out.println("\nAverage number of navigated edges per position " + edgeCount/db.getSize());
    tree = null;
    

    // forward searching
    tree = new KeywordTree(this.queries);
    rootMaxEdge = tree.getRoot().getMaximumEdge().getLabel();
    rootMinEdge = tree.getRoot().getMinimumEdge().getLabel();
    
    msg = String.format("Searching second pass db %d with 1 blind mod: ", db.getSize());
    pm = new ProgressMeter(msg, db.getSize(), System.out);
    edgeCount = 0;
    position = 0;
    
    while (position<db.getSize()) {
      masses.clear();
      
      long dbStart = position;
      while (db.hasMass(position)) {
        masses.add(db.getIntegerMass(position));
        position++;
      }
      
      if (masses.size()>0) {
        Integer[] dbb = masses.toArray(new Integer[masses.size()]);
        
        for (int dbbStart=0; dbbStart<dbb.length; dbbStart++) {
          int firstMass = 0;
          for (int extend=dbbStart; extend<dbb.length; extend++) {
            firstMass += dbb[extend];
            
            // optimization
            if (firstMass < rootMinEdge) continue;
            if (firstMass > rootMaxEdge) break;
            
            int matchIndex = tree.getRoot().search(firstMass);
            if (matchIndex>=0) {
              edgeCount += this.modificationMatch(db, dbStart, dbbStart, dbb, extend+1, tree.getRoot().getEdgeAt(matchIndex), 1, 1, results);
            }
          }
        }
        pm.update(position);
      }
      
      position++;
    }
    System.out.println("\nAverage number of navigated edges per position " + edgeCount/db.getSize());
    tree = null;
  }
  

  
  /**
   * This is the general matching method allowing for 1 mutation.
   * @param db the amino acid sequence 
   * @param dbStart the start position in the db for the first item in dbb
   * @param dbbOffset the number of positions skipped from dbbStart in this search
   * @param cdbb the character db buffer
   * @param dbb the db buffer in integer masses
   * @param dbbIndex the current index to examine in dbb (absolute from the start of dbb)
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param results the results
   * @return the number of edges that we have to navigate
   */
  private int mutationMatch(MassSequence db,
                            long dbStart, 
                            int dbbOffset,
                            Character[] cdbb,
                            Integer[] dbb,
                            int dbbIndex, 
                            Edge edge, 
                            int edgeSubIndex, 
                            int matchedEdgeCount, 
                            ArrayList<MutMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();

    
    // base case because we cannot recurse anymore
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of times we have to call the search function
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      /*
       * we have consumed a complete edge from the compressed tree, and we 
       * arrived to a node, so we need to find a matching edge from the node
       * This is the slow matching situation because we have many branches to match
       */
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps, with no mutations!
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += mutationMatch(db, dbStart, dbbOffset, cdbb, dbb, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, results);
          lowerIndex = matchIndex+1;   // optimization
        }
      }
      
      // for mutation search we have evaluate all possible outgoing edges of this node
      if (matchedEdgeCount>=this.halfLength) {
        
        // different data structures for record keeping
        ArrayList<Mutation> muts = new ArrayList<Mutation>();
        
        for (int edgeIndex=0; edgeIndex<degree; edgeIndex++) {
          Edge targetEdge = sink.getEdgeAt(edgeIndex);
          if (targetEdge==null) continue;
          muts.clear();
          Matching.matchDbWithMutation(dbb, cdbb, dbStart, dbbIndex, targetEdge.getLabel(), muts);
          for (Mutation mutation : muts) {
            edgeCount += exactMatchAfterMut(db, dbStart, dbbOffset, dbb, (int)mutation.getNextStart(), targetEdge, 1, mutation, results);
          }
        }
      }
      
    }
    else {
      /*
       * this case allows for the direct matching without calling the search
       * method of the node because of the unique edge. In other words, we can
       * do a greedy match
       */
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += mutationMatch(db, dbStart, dbbOffset, cdbb, dbb, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, results);
      }
      
      // the mutation search entails searching the only edge against the db
      if (matchedEdgeCount>=this.halfLength) {
        
        // different data structures for record keeping
        ArrayList<Mutation> muts = new ArrayList<Mutation>();
        
        Matching.matchDbWithMutation(dbb, cdbb, dbStart, dbbIndex, currentMass, muts);
        for (Mutation mutation : muts) {
          edgeCount += exactMatchAfterMut(db, dbStart, dbbOffset, dbb, (int)mutation.getNextStart(), edge, edgeSubIndex+1, mutation, results);
        }
      }
      
    }
    
    return edgeCount;
  }
  
  
  
  /**
   * This is the general matching method allowing for 1 mutation for the reverse
   * tree, so we iterate the db sequence in reverse.
   * @param db the amino acid sequence 
   * @param dbStart the start position in the db for the first item in dbb
   * @param dbbOffset the number of positions skipped from dbbStart in this search
   * @param cdbb the character db buffer
   * @param dbb the db buffer in integer masses
   * @param dbbIndex the current index to examine in dbb (absolute from the start of dbb)
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param results the results
   * @return the number of edges that we have to navigate
   */
  private int mutationMatchR(MassSequence db,
                             long dbStart,   
                             int dbbOffset,
                             Character[] cdbb,
                             Integer[] dbb,
                             int dbbIndex,  // the item here has NOT been matched
                             Edge edge, 
                             int edgeSubIndex, 
                             int matchedEdgeCount, 
                             ArrayList<MutMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();
 
    // base case because we don't have anything to recurse
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of times we have to call the search function
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
    	/*
    	 * we have consumed a complete edge from the compressed tree, and we 
       * arrived to a node, so we need to find a matching edge from the node
       * This is the slow matching situation because we have many branches to match
    	 */
    	
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += mutationMatchR(db, dbStart, dbbOffset, cdbb, dbb, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, results);
          lowerIndex = matchIndex+1;   // optimization
        }
      }
      
      // for mutation search we have evaluate all possible outgoing edges of this node
      if (matchedEdgeCount>=this.halfLength) {
      	
        // different data structures for record keeping
      	ArrayList<Mutation> muts = new ArrayList<Mutation>();
        
        for (int edgeIndex=0; edgeIndex<degree; edgeIndex++) {
          Edge targetEdge = sink.getEdgeAt(edgeIndex);
          if (targetEdge==null) continue;
          muts.clear();
          Matching.matchDbWithMutationR(dbb, cdbb, dbStart, dbbIndex, targetEdge.getLabel(), muts);
          for (Mutation mutation : muts) {
            edgeCount += exactMatchAfterMutR(db, dbStart, dbbOffset, dbb, (int)mutation.getNextStart(), targetEdge, 1, mutation, matchedEdgeCount, results);
          }
        }
      }
      
    }
    else {
      // this case allows for the direct matching without calling the search
      // method of the node because of the unique edge. In other words, we can
      // do a greedy match
    	
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += mutationMatchR(db, dbStart, dbbOffset, cdbb, dbb, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, results);
        
      }
      
      // the mutation search entails searching the only edge against the db
      if (matchedEdgeCount>=this.halfLength) {
      	
        // different data structures for record keeping
      	ArrayList<Mutation> muts = new ArrayList<Mutation>();
        Matching.matchDbWithMutationR(dbb, cdbb, dbStart, dbbIndex, currentMass, muts);
        for (Mutation mutation : muts) {
          edgeCount += exactMatchAfterMutR(db, dbStart, dbbOffset, dbb, (int)mutation.getNextStart(), edge, edgeSubIndex+1, mutation, matchedEdgeCount, results);
        }
      }
    }
    
    return edgeCount;
  }
 
  
  
  /**
   * This is the general matching method allowing for 1 modification.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbbOffset the shift from dbb
   * @param dbb the db buffer
   * @param dbbIndex the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param results the results
   * @return the number of edges that we have to navigate
   */
  private int modificationMatch(MassSequence db,
                                long dbStart,
                                int dbbOffset,
                                Integer[] dbb,
                                int dbbIndex, 
                                Edge edge, 
                                int edgeSubIndex, 
                                int matchedEdgeCount, 
                                ArrayList<ModMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();

    
    // base case because we cannot recurse anymore
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of times we have to call the search function
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      // we have consumed a complete edge from the compressed tree, and we 
      // arrived to a node, so we need to find a matching edge from the node
      // This is the slow matching situation because we have many branches to match
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps, with no mutations!
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += modificationMatch(db, dbStart, dbbOffset, dbb, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, results);
          lowerIndex = matchIndex+1;   // optimization
        }
      }
      
      // for mutation search we have evaluate all possible outgoing edges of this node
      if (matchedEdgeCount>=this.halfLength) {
        
        // different data structures for record keeping
        ArrayList<Modification> mods = new ArrayList<Modification>();
        
        for (int edgeIndex=0; edgeIndex<degree; edgeIndex++) {
          Edge targetEdge = sink.getEdgeAt(edgeIndex);
          if (targetEdge==null) continue;
          mods.clear();
          Matching.matchDbWithModification(db, dbStart, dbb, dbbIndex, targetEdge.getLabel(), mods);
          for (Modification mod : mods) {
            edgeCount += exactMatchAfterMod(db, dbStart, dbbOffset, dbb, mod.getNextStart(), targetEdge, 1, matchedEdgeCount+1, mod, results);
          }
        }
      }
      
    }
    else {
      
      // this case allows for the direct matching without calling the search
      // method of the node because of the unique edge. In other words, we can
      // do a greedy match
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += modificationMatch(db, dbStart, dbbOffset, dbb, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, results);
      }
      
      // the mutation search entails searching the only edge against the db
      if (matchedEdgeCount>=this.halfLength) {
        
        // different data structures for record keeping
        ArrayList<Modification> mods = new ArrayList<Modification>();
        
        Matching.matchDbWithModification(db,dbStart, dbb, dbbIndex, currentMass, mods);
        for (Modification mod : mods) {
          edgeCount += exactMatchAfterMod(db, dbStart, dbbOffset, dbb, mod.getNextStart(), edge, edgeSubIndex+1, matchedEdgeCount+1, mod, results);
        }
      }
      
    }
    
    return edgeCount;
  } 
  
  
  
  /**
   * This is the general matching method allowing for 1 modification for the reverse
   * tree, so we iterate the db sequence in reverse.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbbOffset the shift from dbb
   * @param dbb the db buffer
   * @param dbbIndex the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param results the results
   * @return the number of edges that we have to navigate
   */
  private int modificationMatchR(MassSequence db,
                                 long dbStart,
                                 int dbbOffset,
                                 Integer[] dbb,
                                 int dbbIndex, 
                                 Edge edge, 
                                 int edgeSubIndex, 
                                 int matchedEdgeCount, 
                                 ArrayList<ModMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();
 
    
    // base case because we don't have anything to recurse
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of times we have to call the search function
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      /*
       * we have consumed a complete edge from the compressed tree, and we 
       * arrived to a node, so we need to find a matching edge from the node
       * This is the slow matching situation because we have many branches to match
       */
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += modificationMatchR(db, dbStart, dbbOffset, dbb, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, results);
          lowerIndex = matchIndex + 1;
        }
      }
      
      // for mutation search we have evaluate all possible outgoing edges of this node
      if (matchedEdgeCount>=this.halfLength) {
        
        // different data structures for record keeping
        ArrayList<Modification> mods = new ArrayList<Modification>();
        
        for (int edgeIndex=0; edgeIndex<degree; edgeIndex++) {
          Edge targetEdge = sink.getEdgeAt(edgeIndex);
          if (targetEdge==null) continue;
          mods.clear();
          Matching.matchDbWithModificationR(db, dbStart, dbb, dbbIndex, targetEdge.getLabel(), mods);
          for (Modification mod : mods) {
            edgeCount += exactMatchAfterModR(db, dbStart, dbbOffset, dbb, mod.getNextStart(), targetEdge, 1, matchedEdgeCount+1, mod, matchedEdgeCount, results);
          }
        }
      }
      
    }
    else {
      /*
       * this case allows for the direct matching without calling the search
       * method of the node because of the unique edge. In other words, we can
       * do a greedy match
       */
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += modificationMatchR(db, dbStart, dbbOffset, dbb, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, results);
      }
      
      // the mutation search entails searching the only edge against the db
      if (matchedEdgeCount>=this.halfLength) {
        
        // different data structures for record keeping
        ArrayList<Modification> mods = new ArrayList<Modification>();
        Matching.matchDbWithModificationR(db, dbStart, dbb, dbbIndex, currentMass, mods);
        for (Modification mod : mods) {
          edgeCount += exactMatchAfterModR(db, dbStart, dbbOffset, dbb, mod.getNextStart(), edge, edgeSubIndex+1, matchedEdgeCount+1, mod, matchedEdgeCount, results);
        }
      }
    }
    
    return edgeCount;
  }
  
 
  
  /**
   * This is the general matching method that does not allow for any modifications
   * because it has been registered before.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbbOffset the shift from dbb
   * @param dbb the db buffer
   * @param dbbIndex the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param m the mutation object
   * @param results the results
   * @return the number of edges navigated by this search
   */
  private int exactMatchAfterMod(MassSequence db,
                                 long dbStart, 
                                 int dbbOffset,
                                 Integer[] dbb,
                                 int dbbIndex, 
                                 Edge edge, 
                                 int edgeSubIndex, 
                                 int matchedEdgeCount, 
                                 Modification m, 
                                 ArrayList<ModMatchObject> results) {
    
    Node sink = edge.getSink(); 
    int degree = sink.getDegree();  
    
    
    if (edgeSubIndex==edge.size() && sink.getPositionsCount()>0) {
      // there are queries that match with a mutation
      for (int queryIndex : sink.getPositions()) {
        ModMatchObject mo = new ModMatchObject(dbStart+dbbOffset, 
                                               dbStart+dbbIndex+1,
                                               m,
                                               db, 
                                               this.getQueryAt(queryIndex), 
                                               queryIndex);
        results.add(mo);
      }
    }
    
    
    // base case, because cannot keep recursing anymore
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of edges we have to navigate
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
    	// we have consumed a complete edge from the compressed tree, and we 
      // arrived to a node, so we need to find a matching edge from the node
      // This is the slow matching situation because we have many branches to match
    	
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps, with no mutations!
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += exactMatchAfterMod(db, dbStart, dbbOffset, dbb, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, m, results);
          lowerIndex = matchIndex+1;   // optimize the search by resuming at the previous position
        }
      }
    }
    else {
      // this case allows for the direct matching without calling the search
      // method of the node because of the unique edge. In other words, we can
      // do a greedy match
    	
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += exactMatchAfterMod(db, dbStart, dbbOffset, dbb, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, m, results);
      }
    }
    
    return edgeCount;
  }
  
  
  
  /**
   * This is the general matching method that does not allow for any mutations
   * because it has been registered before.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbbOffset the positions shifted from dbb when the current matching was started
   * @param dbb the integer mass buffer
   * @param dbbIndex the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param m the mutation object
   * @param mutationIndex the (edge) index in which the mutation happened (reversed).
   * @param results the results
   * @return the number of edges navigated by this search
   */
  private int exactMatchAfterMut(MassSequence db,
                                 long dbStart, 
                                 int dbbOffset,
                                 Integer[] dbb,
                                 int dbbIndex, 
                                 Edge edge, 
                                 int edgeSubIndex, 
                                 Mutation m, 
                                 ArrayList<MutMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();
    
    
    if (edgeSubIndex==edge.size() && sink.getPositionsCount()>0) {
      // there are queries that match with a mutation
      for (int queryIndex : sink.getPositions()) {
        MutMatchObject smmo = new MutMatchObject(dbStart+dbbOffset, 
                                                 dbStart+dbbIndex,
                                                 m,
                                                 db, 
                                                 this.getQueryAt(queryIndex), 
                                                 queryIndex);
        results.add(smmo);
      }
    }
    
    
    // base case, because cannot keep recursing anymore
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of edges we have to navigate
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      // we have consumed a complete edge from the compressed tree, and we 
      // arrived to a node, so we need to find a matching edge from the node
      // This is the slow matching situation because we have many branches to match
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps, with no mutations!
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += exactMatchAfterMut(db, dbStart, dbbOffset, dbb, i+1, sink.getEdgeAt(matchIndex), 1, m, results);
          lowerIndex = matchIndex+1;   // optimize the search by resuming at the previous position
        }
      }
    }
    else {
      // this case allows for the direct matching without calling the search
      // method of the node because of the unique edge. In other words, we can
      // do a greedy match
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += exactMatchAfterMut(db, dbStart, dbbOffset, dbb, i+1, edge, edgeSubIndex+1, m, results);
      }
    }
    
    return edgeCount;
  }
  
  
  
  /**
   * This is the general matching method that does not allow for any mutations
   * because it has been registered before. This is for reverse matching.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbbOffset the positions shifted from dbb when the current matching was started
   * @param dbb the integer mass buffer
   * @param dbbIndex the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param m the mutation object
   * @param mutationIndex the (edge) index in which the mutation happened (reversed).
   * @param results the results
   * @return the number of edges navigated by this search
   */
  private int exactMatchAfterMutR(MassSequence db,
                                  long dbStart, 
                                  int dbbOffset,
                                  Integer[] dbb,
                                  int dbbIndex,
                                  Edge edge, 
                                  int edgeSubIndex, 
                                  Mutation m,
                                  int mutationIndex,
                                  ArrayList<MutMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();
    
    if (edgeSubIndex==edge.size() && sink.getPositionsCount()>0) {
      // there are queries that match with a mutation
      for (int queryIndex : sink.getPositions()) {
        // we need to check that this have not been added before
        int trueIndex = this.getQueryAt(queryIndex).size()-mutationIndex-1;
        if (!m.isInEdge()) trueIndex++;
        if (trueIndex<this.halfLength) {
          //System.out.printf("start %d, end %d\n", dbStart-dbbIndex, dbStart+1);
          MutMatchObject mo = new MutMatchObject(dbStart-dbbIndex+1, 
                                                 dbStart-dbbOffset+1,
                                                 m,
                                                 db, 
                                                 this.getQueryAt(queryIndex), 
                                                 queryIndex);
          results.add(mo);
        }
      }
    }
    
    
    // base case, because cannot keep recursing anymore
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of edges we navigate
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      /*
       * we have consumed a complete edge from the compressed tree, and we 
       * arrived to a node, so we need to find a matching edge from the node
       * This is the slow matching situation because we have many branches to match
       */
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += exactMatchAfterMutR(db, dbStart, dbbOffset, dbb, i+1, sink.getEdgeAt(matchIndex), 1, m, mutationIndex, results);
          lowerIndex = matchIndex+1; // optimization
        }
      }
    }
    else {
      /*
       * this case allows for the direct matching without calling the search
       * method of the node because of the unique edge. In other words, we can
       * do a greedy match
       */
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += exactMatchAfterMutR(db, dbStart, dbbOffset, dbb, i+1, edge, edgeSubIndex+1, m, mutationIndex, results);
      }
    }
    
    return edgeCount;
  }
  
  
  
  /**
   * This is the general matching method that does not allow for any modifications
   * because it has been registered before. This is for reverse matching.
   * @param db the amino acid sequence 
   * @param dbStart the start position of the database for the current match
   * @param dbbOffset the number of positions shifted on dbb
   * @param dbb the integer mass buffer
   * @param dbbIndex the position to match in the database at this iteration
   * @param edge the edge in the tree to match
   * @param edgeSubIndex the edge subIndex to match because an edge object
   *                     might represent multiple mass edges 
   * @param matchedEdgeCount the number of matched edges so far
   * @param m the modification object
   * @param modIndex the (edge) index in which the modification happened (reversed).
   * @param results the results
   * @return the number of edges navigated by this search
   */
  private int exactMatchAfterModR(MassSequence db,
                                  long dbStart, 
                                  int dbbOffset,
                                  Integer[] dbb,
                                  int dbbIndex, 
                                  Edge edge, 
                                  int edgeSubIndex, 
                                  int matchedEdgeCount, 
                                  Modification m,
                                  int modIndex,
                                  ArrayList<ModMatchObject> results) {
    
    Node sink = edge.getSink();
    int degree = sink.getDegree();
    
   
    if (edgeSubIndex==edge.size() && sink.getPositionsCount()>0) {
      // there are queries that match with a mutation
      for (int queryIndex : sink.getPositions()) {
        // we need to check that this have not been added before
        int trueIndex = this.getQueryAt(queryIndex).size() - modIndex -1;
        if (trueIndex<this.halfLength) {
          ModMatchObject mo = new ModMatchObject(dbStart-dbbIndex+1, 
                                                 dbStart-dbbOffset+1,
                                                 m,
                                                 db, 
                                                 this.getQueryAt(queryIndex), 
                                                 queryIndex);
          results.add(mo);
        }
      }
    }
    
    
    // base case, because cannot keep recursing anymore
    if (edgeSubIndex==edge.size() && degree==0) return 0;
    
    
    int edgeCount = 1;  // counts the number of edges we navigate
    int cumMass = 0;
    if (edgeSubIndex==edge.size()) {
      /*
       * we have consumed a complete edge from the compressed tree, and we 
       * arrived to a node, so we need to find a matching edge from the node
       * This is the slow matching situation because we have many branches to match
       */
      
      // optimization stuff
      int lower = sink.getMinimumEdge().getLabel();
      int upper = sink.getMaximumEdge().getLabel();
    
      // branch many possible gaps
      int lowerIndex = 0;
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];

        // optimization
        if (cumMass < lower) continue;
        if (cumMass > upper) break;
        
        int matchIndex = sink.search(cumMass, lowerIndex, degree);
        if (matchIndex>=0) {
          // recurse
          edgeCount += exactMatchAfterModR(db, dbStart, dbbOffset, dbb, i+1, sink.getEdgeAt(matchIndex), 1, matchedEdgeCount+1, m, modIndex, results);
          lowerIndex = matchIndex + 1; // optimization
        }
      }
    }
    else {
      /*
       * this case allows for the direct matching without calling the search
       * method of the node because of the unique edge. In other words, we can
       * do a greedy match
       */
      
      int currentMass = edge.getLabelAt(edgeSubIndex);
      for (int i=dbbIndex; i<dbb.length; i++) {
        cumMass += dbb[i];
        
        if (cumMass < currentMass) continue;
        // we have reached an impossible path, by over shooting
        if (cumMass > currentMass) break;
        
        // There is a match
        edgeCount += exactMatchAfterModR(db, dbStart, dbbOffset, dbb, i+1, edge, edgeSubIndex+1, matchedEdgeCount+1, m, modIndex, results);
      }
    }
    
    return edgeCount;
  }

  
} // end of class
