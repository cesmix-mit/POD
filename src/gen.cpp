#include "Halide.h"
#include <stdio.h>
 
using namespace Halide;

//Func buildRbfFunc(std::string call, Func rij, Func scalefunc, Expr rin, Expr rmax, int pdegree, int K, )

std::tuple<Func, Func> buildNeighPairs(std::string call, Func pairnumsum, Func pairlist,
				       Expr NPairs, Expr NAtoms, Expr NMax, Expr dim, Expr NTypes,
				       Func atomtype, Func alist, Func atompos
		     ){


  Func outputs(call + "_bnp");
  Func vectors(call + "_bnpd");

  Var ii("ii_" + call), c("c_" + call), j("j_" + call);
  //initialize ouputs
  outputs(j, c) = mux(c, {-1, -1, -1, -1});//ai, aj, ti, tj
  outputs.bound(c, 0, 4);
  outputs.bound(j, 0, NPairs);

  vectors(j, c) = Expr((double) 0.0);
  vectors.bound(c, 0, dim);
  vectors.bound(j, 0, NPairs);


  RDom r(0, NAtoms, 0, NMax);
  r.where(r.y < pairnumsum(r.x + 1) && r.y >= pairnumsum(r.x)); //iteration

  Expr jacc = clamp(pairlist(r.y), 0, NPairs);
  Expr att = clamp(alist(jacc), 0, NAtoms); 
  Expr att_tt = atomtype(att); 
  outputs(r.x, c) = mux(c, {r.x, att, atomtype(r.x), att_tt});

  vectors(r.x, c) = atompos(r.y, c) - atompos(r.x, c);

  return std::make_tuple(outputs, vectors);

  
}

class pod1 : public Halide::Generator<pod1> {
public:
  //Func pairnumsum, Func pairlist,
  //				       Expr NPairs, Expr NAtoms, Expr NMax, Expr dim,
  //				       Func atomtype, Func alist, Func atompos
  
  Input<Buffer<int>> pairlist{"pairlist", 1};
  Input<Buffer<int>> pairnumsum{"pairnumsum", 1};
  Input<Buffer<int>> atomtype{"atomptype", 1};
  Input<Buffer<int>> alist{"alist", 1};
  Input<Buffer<double>> atompos{"atompos", 2};

  Output<Buffer<double>> rij{"rij", 2};
  Output<Buffer<int>> meta{"meta", 2};

  Pipeline pipeline;

  GeneratorParam<int> NMax{"NMax", 100};
  GeneratorParam<int> NTypes{"NTypes", 3};
  //  GeneratorParam<int> M{"M", 50};
  //  GeneratorParam<int> I{"I", 5000};
  


  void generate (){

    Expr NPairs = pairlist.dim(0).max();
    Expr NAtoms = atomtype.dim(0).max();
    alist.dim(0).set_bounds(0, NAtoms);
    atompos.dim(0).set_bounds(0, NAtoms);
    atompos.dim(1).set_bounds(0, 3);
    Func np1_data, np1_vecs;
    std::tie(np1_data, np1_vecs) = buildNeighPairs("np1", pairnumsum, pairlist,
						   NPairs, NAtoms, NMax, 3, NTypes,
						     atomtype, alist, atompos);

    Var i("i"), j("j");

    rij(i, j)= np1_vecs(i, j);
    meta(i, j) = np1_data(i, j);

  }
    
    

    



};
  
  

HALIDE_REGISTER_GENERATOR(pod1, pod1);
// uquad1(qy, qx, y, x) = qwx(qx) * qwy(qy) * sum(rpx, basisx(qx, rpx.x) * sum(rpy, basisyp(qy, rpy.x) * dist(rpy.x, rpx.x, y, x)));
// uquad2(qy, qx, y, x) = qwx(qx) * qwy(qy) * sum(rpx, basisxp(qx, rpx.x) * sum(rpy, basisy(qy, rpy.x) * dist(rpy.x, rpx.x, y, x)));
// iquad(yp, xp, y, x) = sum(r, uquad1(r.y, r.x, y, x) * basisx(r.x, xp) * basisyp(r.y, yp) + uquad2(r.y, r.x, y, x) * basisxp(r.x, xp) * basisy(r.y, yp));
// compute(yp, xp, y, x) = iquad(yp, xp, y, x);
// qwx(qx) * qwy(qy) * sum(rpx, basisx(qx, rpx.x) * sum(rpy, basisyp(qy, rpy.x) * dist(rpy.x, rpx.x, y, x)));

