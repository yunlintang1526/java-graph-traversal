/*
 * NAME: Yunlin Tang
 * PID: a14664383
 */

import java.util.*;


/**
 * the graph class which construct a graph includes vertex and edges
 *
 * @author Yunlin Tang
 * @since ${12.3}
 */
public class Graph {

    private HashMap graph; // the main map to store the vertex
    private List<Edge> path; // path that connect two vertex
    private List<Edge> allEdges; // the list store all the edges on the graph
    private final static int SQUARE = 2; // square of the number


    /**
     * Constructor for Graph
     */
    public Graph() {
        // initialize all the useful data structures to store information
        graph = new HashMap();
        path = new ArrayList<>();
        allEdges = new ArrayList<>();
    }

    /**
     * Adds a vertex to the graph. Throws IllegalArgumentException if given vertex
     * already exist in the graph.
     *
     * @param v vertex to be added to the graph
     * @throws IllegalArgumentException if two vertices with the same name are added.
     */
    @SuppressWarnings("unchecked")
    public void addVertex(Vertex v) throws IllegalArgumentException {
        // if the graph already contains the vertex, throw exception
        if(graph.containsKey(v.hashCode())){
            throw new IllegalArgumentException();
        }
        //add the new vertex to the map
        graph.put(v.hashCode(), v);
    }

    /**
     * Gets a collection of all the vertices in the graph
     *
     * @return collection of all the vertices in the graph
     */
    @SuppressWarnings("unchecked")
    public Collection<Vertex> getVertices() {
        return graph.values();
    }

    /**
     * Gets the vertex object with the given name
     *
     * @param name name of the vertex object requested
     * @return vertex object associated with the name
     */
    public Vertex getVertex(String name) {
        // loop over the hashmap's entryset
        Iterator iterator = graph.values().iterator();
        while(iterator.hasNext()){
            Vertex v = (Vertex)iterator.next();
            // check if the current vertex's name equal to the input name
            if(v.getName().equals(name)){
                return v;
            }
        }
        return null;
    }

    /**
     * Adds a directed edge from vertex u to vertex v, Throws
     * IllegalArgumentException if one of the vertex does not exist
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight weight of the edge between vertex u and v
     * @throws IllegalArgumentException if one of the vertex does not exist
     */
    public void addEdge(String nameU, String nameV, Double weight)
            throws IllegalArgumentException {

        // get the corresponding vertex from the string name
        Vertex u = getVertex(nameU);
        Vertex v = getVertex(nameV);

        // if one of them is null, throws exception
        if(u == null || v == null){
            throw new IllegalArgumentException();
        }

        // set the u's children to v
        u.setNext(v,weight);

        // add the this new edge to the edge list
        allEdges.add(new Edge(u,v,weight));
    }

    /**
     * Adds an undirected edge between vertex u and vertex v by adding a directed
     * edge from u to v, then a directed edge from v to u
     *
     * @param nameU name of vertex u
     * @param nameV name of vertex v
     * @param weight  weight of the edge between vertex u and v
     */
    public void addUndirectedEdge(String nameU, String nameV, double weight) {

        // call two times of addEdge to add undirected edge
        addEdge(nameU, nameV, weight);
        addEdge(nameV, nameU, weight);


    }

    /**
     * getter to get all the edges in the graph
     * @return the list stored all the edges
     */
    public List<Edge> getAllEdges(){
        return allEdges;
    }


    /**
     * Computes the euclidean distance between two points as described by their
     * coordinates
     *
     * @param ux (double) x coordinate of point u
     * @param uy (double) y coordinate of point u
     * @param vx (double) x coordinate of point v
     * @param vy (double) y coordinate of point v
     * @return (double) distance between the two points
     */
    public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {

        return Math.sqrt(Math.pow(vx-ux,SQUARE)+Math.pow(vy-uy, SQUARE));
    }

    /**
     * Calculates the euclidean distance for all edges in the map using the
     * computeEuclideanCost method.
     */
    @SuppressWarnings("unchecked")
    public void computeAllEuclideanDistances() {
        // loop over all the vertex
        Iterator iterator = graph.values().iterator();
        while(iterator.hasNext()){
            Vertex v = (Vertex)iterator.next();
            // loop over all the children of the current node
            ArrayList<Edge> edge = v.getNext();
            for(Edge e: edge){
                Vertex source = e.getSource();
                Vertex target = e.getTarget();
                // calculate the distance
                e.setDistance(computeEuclideanDistance(source.getX(),source.getY(),
                        target.getX(), target.getY()));
            }
        }

    }

    /**
     * Helper method to reset all the vertices before doing graph traversal algorithms
     */
    private void resetAllVertices() {
        // reset the previous variable of the current node
        for(Object v:graph.values()){
            ((Vertex)v).setPrev(null);
        }
        // reset the path
        path = new ArrayList<>();

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using DFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void DFS(String s, String t) {

        resetAllVertices();
        // initialize all the data structure that relevant to this method
        Stack<Vertex> stack = new Stack<>();
        ArrayList<Vertex> visited = new ArrayList<>();
        Vertex start = getVertex(s);
        Vertex end = getVertex(t);
        stack.push(start);

        // if two nodes are equal, just return
        if(start.equals(end)){
            path.add(new Edge(start, end, 0));
            return;
        }

        // loop on the stack as long as it is not empty
        while(!stack.isEmpty()){
            // get the node at the top of the stack
            Vertex currV = stack.pop();
            // make it visited
            visited.add(currV);
            // loop over every child of this current node
            for(Edge e: currV.getNext()){
                // if the child is not visited yet
                if(!visited.contains(e.getTarget())) {
                    // if the child equals to the end node, return
                    if (e.getTarget().equals(end)) {
                        end.setPrev(currV);
                        stack.clear();
                        break;
                    }
                    // otherwise set the previous of child to the current node
                    // and push it to the stack
                    e.getTarget().setPrev(currV);
                    stack.push(e.getTarget());
                }
            }
        }

        // trace back the path which starts from the end node
        Stack<Edge> temp = new Stack<>();
        Vertex currV = end;
        while(!currV.equals(start)){
            // search on the edge list
            for(Edge e: allEdges){
                if(e.getTarget().equals(currV) &&
                        e.getSource().equals(currV.getPrev())){
                    temp.push(e);
                }
            }
            currV = currV.getPrev();
        }

        // add the edge to the path
        while(!temp.isEmpty()){
            path.add(temp.pop());
        }

    }

    /**
     * Find the path from vertex with name s to vertex with name t, using BFS
     *
     * @param s the name of the starting vertex
     * @param t the name of the targeting vertex
     */
    public void BFS(String s, String t) {

        resetAllVertices();
        // initialize all the data structure which needed in this method
        LinkedList<Vertex> frontier = new LinkedList<>();
        ArrayList<Vertex> visited = new ArrayList<>();
        Vertex start = getVertex(s);
        Vertex end = getVertex(t);
        // add the root to the frontier list
        frontier.add(start);

        // if the end node and start node equals, just return
        if(start.equals(end)){
            path.add(new Edge(start, end, 0));
            return;
        }

        // mark the root as visited
        visited.add(start);
        // loop on the frontier list as long as it is not empty
        while(!frontier.isEmpty()){
            // get and remove the node at the top of the list
            Vertex currV = frontier.poll();
            // get the children of current vertex
            for(Edge e: currV.getNext()){
                // if the child is not yet visited, continue
                if(!visited.contains(e.getTarget())){
                    // mark the child is visited
                    visited.add(e.getTarget());
                    if(e.getTarget().equals(end)){
                        // if the child is end vertex, set the previous
                        // then stop the iteration
                        end.setPrev(currV);
                        frontier.clear();
                        break;
                    }
                    // otherwise set the previous, and push it to the list
                    e.getTarget().setPrev(currV);
                    frontier.add(e.getTarget());
                }
            }
        }

        // trace back the path, which starts from the end vertex
        Stack<Edge> temp = new Stack<>();
        Vertex currV = end;
        // loop on the path
        while(!currV.equals(start)){
            for(Edge e: allEdges){
                if(e.getTarget().equals(currV) &&
                        e.getSource().equals(currV.getPrev())){
                    temp.push(e);
                }
            }
            currV = currV.getPrev();
        }

        // add the vertex from the stack to the path
        while(!temp.isEmpty()){
            path.add(temp.pop());
        }

    }

    /**
     * Helper class for Dijkstra and A*, used in priority queue
     */
    private class CostVertex implements Comparable<CostVertex> {
        // all the relevant variables to the class
        double cost;
        Vertex vertex;
        CostVertex prev; // previous vertex from the path

        /**
         * the constructor of CostVertex
         * @param cost the distance from the current node to the source
         * @param vertex the vertex which this CostVertex includes
         */
        public CostVertex(double cost, Vertex vertex) {
            this.cost = cost;
            this.vertex = vertex;
            prev = null;
        }

        /**
         * compare method between this CostVertex to other CostVertex
         * @param o other CostVertex
         * @return the integer indicate which is greater
         */
        public int compareTo(CostVertex o) {
            return Double.compare(cost, o.cost);
        }

        /**
         * setter method for previous variable
         * @param p the previous vertex from the path
         */
        public void setPrev(CostVertex p){
            prev = p;
        }

        /**
         * the getter method for previous variables
         * @return the previous vertex
         */
        public CostVertex getPrev(){
            return prev;
        }
    }

    /**
     * Find the shortest path from vertex with name s to vertex with name t, using Dijkstra
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    @SuppressWarnings("unchecked")
    public void Dijkstra(String s, String t) {
        resetAllVertices();
        // initialize all the relevant data structure
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        Collection<Vertex> allVertex = graph.values();
        ArrayList<CostVertex> allCost = new ArrayList<>();
        // get the source and target vertex
        Vertex start = getVertex(s);
        Vertex end = getVertex(t);

        // store all the CostVertex to the allCost list
        for(Vertex v: allVertex){
            if(v.equals(start)){
                // assign root distance to 0
                allCost.add(new CostVertex(0, start));
            } else {
                // assign other vertex distance to max integer value
                allCost.add(new CostVertex(Integer.MAX_VALUE, v));
            }
        }
        // add all the vertex to the queue
        queue.addAll(allCost);

        // loop on the queue as long as it is not empty
        while(!queue.isEmpty()){
            // choose the vertex with the minimum distance from the queue
            CostVertex currV = queue.poll();

            // if the current vertex is the end vertex, return
            if(currV.vertex.equals(end)){
                queue.clear();
                break;
            }

            // loop on every child of current vertex
            for(Edge e:currV.vertex.getNext()){
                // determine if the child is visited or not yet
                boolean notVisited = false;
                for(CostVertex c: queue){
                    if(c.vertex.equals(e.getTarget())){
                        notVisited = true;
                        break;
                    }
                }

                // if the child is not visited
                if(notVisited) {
                    // get the distance from the current cost and the distance
                    double tempDist = currV.cost + e.getDistance();

                    // get the corresponding costVertex
                    for (CostVertex v : allCost) {
                        if (v.vertex.equals(e.getTarget())){
                            // if the temp distance is smaller than the current
                            // distance(cost), assign it
                            if (tempDist < v.cost) {
                                v.cost = tempDist;
                                v.setPrev(currV);
                                // replace the edited CostVertex
                                for(CostVertex c: queue){
                                    if(c.vertex.equals(v.vertex)){
                                        queue.remove(c);
                                        queue.add(v);
                                        break;
                                    }
                                }
                            }
                        }
                    }


                }
            }
        }

        // GET the current vertex from the CostVertex
        CostVertex currV = null;
        for(CostVertex v: allCost){
            if(v.vertex.equals(end)){
                currV = v;
            }
        }
        // trace back the path
        Stack<Edge> temp = new Stack<>();
        while(!currV.vertex.equals(start)){
            for(Edge e: allEdges){
                if(e.getSource().equals(currV.getPrev().vertex) &&
                    e.getTarget().equals(currV.vertex)){
                    temp.add(e);
                }
            }
            currV = currV.getPrev();
        }

        while(!temp.isEmpty()){
            path.add(temp.pop());
        }

    }

    /**
     * Helper method to calculate the h value in A*
     *
     * @param cur the current vertex being explored
     * @param goal the goal vertex to reach
     * @return the h value of cur and goal vertices
     */
    private double hValue(String cur, String goal) {

        Vertex currV = getVertex(cur);
        Vertex goalV = getVertex(goal);
        // calculate the h-value
        return computeEuclideanDistance(currV.getX(), currV.getY(),
                    goalV.getX(), goalV.getY());
    }

    /**
     * Find the path from vertex with name s to vertex with name t, using A*
     *
     * @param s the name of starting vertex
     * @param t the name of targeting vertex
     */
    @SuppressWarnings("unchecked")
    public void AStar(String s, String t) {

        resetAllVertices();
        // initialize all the data structure
        PriorityQueue<CostVertex> queue = new PriorityQueue<>();
        Collection<Vertex> allVertex = graph.values();
        ArrayList<CostVertex> allCost = new ArrayList<>();

        Vertex start = getVertex(s);
        Vertex end = getVertex(t);

        // add the costVertex to the list corresponding to the vertex on the graph
        for(Vertex v: allVertex){
            if(v.equals(start)){
                allCost.add(new CostVertex(0+hValue(s,t), start));
            } else {
                allCost.add(new CostVertex(Integer.MAX_VALUE, v));
            }
        }
        queue.addAll(allCost);

        // loop on the queue as long as it is not empty
        while(!queue.isEmpty()){
            CostVertex currV = queue.poll();

            if(currV.vertex.equals(end)){
                queue.clear();
                break;
            }
            // loop on all the edge which are the children of the current vertex
            for(Edge e:currV.vertex.getNext()){
                // determine if the child is visited or not
                boolean notVisited = false;
                for(CostVertex c: queue){
                    if(c.vertex.equals(e.getTarget())){
                        notVisited = true;
                        break;
                    }
                }
                // if the child is not visited yet
                if(notVisited) {
                    // get the temp distance by calculating the h-value
                    double tempDist = currV.cost + e.getDistance() +
                            hValue(e.getTarget().getName(), t) -
                            hValue(e.getSource().getName(),t);
                    // update the vertex
                    for (CostVertex v : allCost) {
                        if (v.vertex.equals(e.getTarget())){
                            if (tempDist < v.cost) {
                                v.cost = tempDist;
                                v.setPrev(currV);
                                // replace the vertex in the list
                                for(CostVertex c: queue){
                                    if(c.vertex.equals(v.vertex)){
                                        queue.remove(c);
                                        queue.add(v);
                                        break;
                                    }
                                }
                            }
                        }
                    }


                }
            }
        }

        CostVertex currV = null;
        for(CostVertex v: allCost){
            if(v.vertex.equals(end)){
                currV = v;
            }
        }

        Stack<Edge> temp = new Stack<>();
        while(!currV.vertex.equals(start)){
            for(Edge e: allEdges){
                if(e.getSource().equals(currV.getPrev().vertex) &&
                        e.getTarget().equals(currV.vertex)){
                    temp.push(e);
                }
            }
            currV = currV.getPrev();
        }
        while(!temp.isEmpty()){
            path.add(temp.pop());
        }

    }

    /**
     * Returns a list of edges for a path from city s to city t.
     *
     * @param s starting city name
     * @param t ending city name
     * @return list of edges from s to t
     */
    public List<Edge> getPath(String s, String t) {
        return path;
    }

}
