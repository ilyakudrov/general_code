#include "link.h"
#include "monopoles.h"

loop::loop (const link1<matrix> n1) {node=link1<matrix>(n1); for(int i=0; i<=6; i++) link[i]=NULL;};
loop::loop() { node=link1<matrix>(); for(int i=0; i<=6; i++) link[i]=NULL; };
loop::loop( const loop& l) { node=link1<matrix>(l.node);  for(int i=0; i<=6; i++) link[i]=l.link[i]; };

int loop::get_dir(int i) {
      loop* tmp=link[i-1];
      if(tmp==NULL) return 0;

      int delta=(tmp->node.coordinate[0])-node.coordinate[0];
      if(delta==(1-x_size)) return 1; if(delta==(x_size-1)) return -1;
      if(delta!=0) return 1*delta;

      delta=(tmp->node.coordinate[1])-node.coordinate[1];
      if(delta==(1-y_size)) return 2; if(delta==(y_size-1)) return -2;
      if(delta!=0) return delta*2;

      delta=(tmp->node.coordinate[2])-node.coordinate[2];
      if(delta==(1-z_size)) return 3; if(delta==(z_size-1)) return -3;
      if(delta!=0) return delta*3;

      delta=(tmp->node.coordinate[3])-node.coordinate[3];
      if(delta==(1-t_size)) return 4; if(delta==(t_size-1)) return -4;
      if(delta!=0) return delta*4;
      return 0;
}

void calculate_current(data_matrix& conf, FLOAT* J){
  link1<matrix> link(x_size, y_size, z_size, t_size);
	int a;
	for(int t=1;t<=t_size;t++) {
    for(int z=1;z<=z_size;z++) {
      for(int y=1;y<=y_size;y++) {
        for(int x=1;x<=x_size;x++) {
          link.go(x, y, z, t);
          for(int i = 1;i <= 4;i++) {
            link.move_dir(i);
            a = link.get_place();
            link.move(i, 1);
            J[a]=link.get_current(conf);
            link.move(i, -1);
          }
				}
			}
		}
	}
}

int test(link1<matrix> l, FLOAT* J){
  for(int i=1; i<=4; i++){
    if( (J[l.get_place1()*4+i-1]>0.3)||(J[l.get_place1()*4+i-1]<-0.3) ) return i;
  }
  for(int i=1; i<=4; i++){
    l.move(i,-1);
    if( (J[l.get_place1()*4+i-1]>0.3)||(J[l.get_place1()*4+i-1]<-0.3) ) return -i;
    l.move(i,1);
  }
  return 0;
}

void find_cluster(loop* ll, FLOAT* J) {
  int i=0;
  link1<matrix> tmp;
  int dir;
  do {
    dir=test(ll->node, J);
    tmp=ll->node;
    if (dir>0)  { J[(ll->node.get_place1()*4)+dir-1]=0.; tmp.move(dir,1); };
    if (dir<0)  {link1<matrix> tmp1(ll->node); tmp1.move(-dir,-1);
      J[(tmp1.get_place1()*4)-(dir)-1]=0.; tmp.move(-dir,-1);};

    loop* l1=new loop(tmp);
    ll->link[i]=l1;
    i++;
    dir=test(tmp, J);
    if(dir==0) continue;
    else find_cluster(l1, J);
  } while(test(ll->node, J)!=0);
  return ;
};

void calculate_clusters(vector<loop*>& LL, FLOAT* J){
  int dir1;
  for(int it=1;it<=t_size; it++ ) {
    for(int iz=1;iz<=z_size; iz++ ) {
      for(int iy=1;iy<=y_size; iy++ ) {
        for(int ix=1;ix<=x_size; ix++ ) {
          link1<matrix> tmp(ix,iy,iz,it, x_size, y_size, z_size, t_size);
          dir1=test(tmp, J);
          if(dir1==0) continue;
          LL.push_back(new loop(tmp));
          find_cluster( LL[LL.size()-1], J );
        }
      }
    }
  }
}
