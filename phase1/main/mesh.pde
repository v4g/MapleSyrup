// TRIANGLE MESH
class MESH {
    // VERTICES
    int nv=0, maxnv = 1000;  
    pt[] G = new pt [maxnv];                        
    pt[] cc = new pt [maxnv];
    // TRIANGLES 
    int nt = 0, maxnt = maxnv*2;                           
    boolean[] isInterior = new boolean[maxnv];                                      
    // CORNERS 
    int c=0;    // current corner                                                              
    int nc = 0; 
    int[] V = new int [3*maxnt];   
    int[] O = new int [3*maxnt];  
    
    // current corner that can be edited with keys
  MESH() {for (int i=0; i<maxnv; i++) G[i]=new pt();};
  void reset() {nv=0; nt=0; nc=0;}                                                  // removes all vertices and triangles
  void loadVertices(pt[] P, int n) {nv=0; for (int i=0; i<n; i++) addVertex(P[i]);}
  void writeVerticesTo(pts P) {for (int i=0; i<nv; i++) P.G[i].setTo(G[i]);}
  void addVertex(pt P) { G[nv++].setTo(P); }                                             // adds a vertex to vertex table G
  void addTriangle(int i, int j, int k) {V[nc++]=i; V[nc++]=j; V[nc++]=k; nt=nc/3; }     // adds triangle (i,j,k) to V table
  void addCircumcenter(int i, pt p) {cc[i]=p;}     // adds triangle (i,j,k) to V table
  void printTriangle(int i){println(V[i%nc],V[(i+1)%nc],V[(i+2)%nc]);}
  // CORNER OPERATORS
  int t (int c) {int r=int(c/3); return(r);}                   // triangle of corner c
  int n (int c) {int r=3*int(c/3)+(c+1)%3; return(r);}         // next corner
  int p (int c) {int r=3*int(c/3)+(c+2)%3; return(r);}         // previous corner
  pt g (int c) {return G[V[c]];}                             // shortcut to get the point where the vertex v(c) of corner c is located
  pt cc (int c) {return cc[t(c)];}
  boolean nb(int c) {return(O[c]!=c);};  // not a border corner
  boolean bord(int c) {return(O[c]==c);};  // not a border corner
  int coft (int c) {return 3*c;} //find first corner of triangle
  pt cg(int c) {return P(0.6,g(c),0.2,g(p(c)),0.2,g(n(c)));}   // computes offset location of point at corner c

  // CORNER ACTIONS CURRENT CORNER c
  void next() {c=n(c);}
  void previous() {c=p(c);}
  void opposite() {c=o(c);}
  void left() {c=l(c);}
  void right() {c=r(c);}
  void swing() {c=s(c);} 
  void unswing() {c=u(c);} 
  void printCorner() {println("c = "+c);}
  
  

  // DISPLAY
  void showCurrentCorner(float r) { if(bord(c)) fill(red); else fill(dgreen); show(cg(c),r); };   // renders corner c as small ball
  void showEdge(int c) {beam( g(p(c)),g(n(c)),rt ); };  // draws edge of t(c) opposite to corner c
  void showVertices(float r) // shows all vertices green inside, red outside
    {
    for (int v=0; v<nv; v++) 
      {
      if(isInterior[v]) fill(green); else fill(red);
      show(G[v],r);
      }
    }                          
  void showInteriorVertices(float r) {for (int v=0; v<nv; v++) if(isInterior[v]) show(G[v],r); }                          // shows all vertices as dots
  void showTriangles() { for (int c=0; c<nc; c+=3) /*if(triangleN%nt == c/3)*/show(g(c), g(c+1), g(c+2)); }         // draws all triangles (edges, or filled)
  void showEdges() {for (int i=0; i<nc; i++) showEdge(i); };         // draws all edges of mesh twice

  void triangulate()      // performs Delaunay triangulation using a quartic algorithm
  {
   c=0;                   // to reset current corner
   for (int i = 0; i < nv; i++)
   {
     for (int j = i+1; j < nv; j++)
     {
       for (int k = j+1; k < nv ; k++)
       {
         pt cc = CircumCenter(G[i],G[j],G[k]);
         float r = d(G[i],cc);
         int l = 0;
         for (l = 0; l < nv ; l++)
           if (l != i & l != j & l != k && d(G[l],cc) <= r) break;
         if(l == nv)
         { 
           if(cw(V(0,0,1),V(G[j],G[k]),V(G[k],G[i])))
             addTriangle(i,j,k);
           else
             addTriangle(k,j,i);
           addCircumcenter(nt-1,cc);
         }
       }
     }
   }
 }  

   
  void computeO() // **02 implement it 
  { 
    for (int i = 0 ; i < nc ; i++)  O[i] = i;
    for (int i = 0 ; i < nc ; i++)
    {
      for (int j = i+1; j < nc; j++)
      {
        if(( v(n(i)) == v(p(j)) ) && ( v(p(i)) == v(n(j)) )){  
          O[i] = j;O[j] = i;
          }  
      }
    }
  } 
    
  void showBorderEdges()  // draws all border edges of mesh
  {
    for (int i = 0 ; i < nc; i++) if(bord(i)) showEdge(i);
  }         

  void showNonBorderEdges() // draws all non-border edges of mesh
  {
    for (int i = 0 ; i < nc; i++) if(nb(i)) showEdge(i);
  }        
    
  void classifyVertices() 
  {
    for (int i = 0; i < nv; i++) isInterior[i] = true;
    for (int i = 0; i < nc; i++){
      if(bord(i)){ //this is exterior
        isInterior[v(n(i))] = false; isInterior[v(p(i))] = false;
      }
    }
  }  
    
  void smoothenInterior() { // even interior vertiex locations
    pt[] Gn = new pt[nv];
    boolean[] visited = new boolean[nv];
    for (int i=0; i<nv; i++) visited[i] = false;
    for (int i=0; i<nc; i++){
      if(!visited[v(i)] && isInterior[v(i)]){
        int c = i;float count=0;
        pt sum = P();
        do{
          sum = A(sum,P(triArea(c),triCenter(c))); count += triArea(c);
          c = s(c);
        }while(c != i);
        Gn[v(i)] = P(1./count,sum);
      }
    }
    for (int v=0; v<nv; v++) if(isInterior[v]) G[v].translateTowards(.1,Gn[v]);
    }


   // **05 implement corner operators in Mesh
  int v (int c) {return V[c];}                                // vertex of c
  int o (int c) {return O[c];}                                // opposite corner
  int l (int c) {return o(n(c));}                             // left
  int s (int c) {return n(o(n(c)));}                             // left
  int u (int c) {return p(o(p(c)));}                             // left
  int r (int c) {return o(p(c));}                             // right

  void showVoronoiEdges() // draws Voronoi edges on the boundary of Voroni cells of interior vertices
  {
    for(int i=0; i<nt; i++){
      int c = coft(i); //corner of this triangle
      //println("triangle ",i," corner ",c," opposite ",o(c));
      show(cc(c),cc(o(c)));
      show(cc(c),cc(o(n(c))));
      show(cc(c),cc(o(p(c))));
    }
  }               

  void drawBezier(pt a, pt b, pt c)
  {
    float dt = 1./20;
    beginShape(LINE);
    for (float t = 0.0; t < 1.0f; t+=dt)show(Bezier(a,b,c,t),Bezier(a,b,c,t+dt));
    endShape();
  }
  void showOpposites() // draws Voronoi edges on the boundary of Voroni cells of interior vertices
  {
    for (int i=0;i<nc;i++)
      if(nb(i)) drawBezier(g(i),P(g(p(i)),g(n(i))),g(o(i)));
  }               

  void showArcs() // draws arcs of quadratic B-spline of Voronoi boundary loops of interior vertices
  {
    for(int i=0; i<nt; i++){
      int c = coft(i); //corner of this triangle
      for (int j = 0; j < 3; j++){
        drawBezier(P(cc(o(c)),cc(c)),cc(c),P(cc(o(n(c))),cc(c)));
        c = n(c);  
      }
    }
  }               // draws arcs in triangles

 
  pt triCenter(int c) {return P(g(c),g(n(c)),g(p(c))); }  // returns center of mass of triangle of corner c
  pt triCircumcenter(int c) {return CircumCenter(g(c),g(n(c)),g(p(c))); }  // returns circumcenter of triangle of corner c
  float triArea(int c) { return 0.5*det3(V(g(c),g(n(c))),V(g(c),g(p(c)))); }


  } // end of MESH
