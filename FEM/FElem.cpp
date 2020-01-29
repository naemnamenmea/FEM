////------------------------------------------------------------------------------
//#include "FElem.h"
//#include "Analysis.h"
//#include "Material.h"
//#ifdef STRONGCHECK
// #include "ThrowMessage.h"
//#endif
////------------------------------------------------------------------------------
//#include <algorithm>
////------------------------------------------------------------------------------
//using namespace std;
////------------------------------------------------------------------------------
////--------------------------tFinitElement---------------------------------------
////------------------------------------------------------------------------------
//static const cardinal_t number_of_FE_kinds = 3;
//static const char*         initNames[]    = {"Rod2Const", "MembraneQuad4Const", "Parallelepiped8"};//Add new kind names here
//static tFinitElement*(*initFunctions[])() = {tIsoRod2ConstSec::NewFE,tIsoQuad4ConsThick::NewFE,tIsoParallelepiped8::NewFE};//Add new pointers to creaters here
//
//const tNamesToFuns<tFinitElement> tFinitElement::Factory(initNames,initFunctions,number_of_FE_kinds);
//
//const string&    tIsoRod2ConstSec::KindName = tFinitElement::Factory.FindName(tIsoRod2ConstSec::NewFE);
//const string&  tIsoQuad4ConsThick::KindName = tFinitElement::Factory.FindName(tIsoQuad4ConsThick::NewFE);
//const string& tIsoParallelepiped8::KindName = tFinitElement::Factory.FindName(tIsoParallelepiped8::NewFE);
////------------------------------------------------------------------------------
//void tFinitElement::DefineNextNode (const tNode& node_)
//{
//#ifdef STRONGCHECK
// std::vector<const tNode*>::iterator ppnode = std::find(pNodes.begin(),pNodes.end(),(const tNode*)NULL);
// Assert(ppnode < pNodes.end(), "Attempt to link redundant node to FE");
//#endif
// *std::find(pNodes.begin(),pNodes.end(),(const tNode*)NULL) = &node_;
//}
////------------------------------------------------------------------------------
//const SymmetricTensor4s& tFinitElement::MaterialElasTensor (const Tensor1s& lc_) const
//{
// return pMaterial->ElasTensor(TensorOfRotationToLocal(lc_),StrainLevel);
//}
////------------------------------------------------------------------------------
//Tensor1s& tFinitElement::GlobalCoord (const Tensor1s& lc_, Tensor1s& result_) const
//{
// result_.Assign0();   Tensor1s radiusvector;   const tNode* pnode;
// for (cardinal_t i=1; i<=HowManyCoordApproxNodes(); ++i)
//   {
//    pnode = &CoordApproxNode(i);
//    radiusvector = pnode->Coord();
//    result_ += (radiusvector *= CoordShapeFun(*pnode,lc_));
//   }
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor2s& tFinitElement::JacobyMatrix_l2g (const Tensor1s& lc_, Tensor2s& result_) const
//{
// result_.Assign0();   Tensor1s tmp1;   Tensor2s tmp2;   const tNode* pnode;
// for (cardinal_t i=1; i<=HowManyCoordApproxNodes(); ++i)
//   {
//    pnode = &CoordApproxNode(i);
//    result_ += CoordShapeGrad(*pnode,lc_,tmp1).DirectProduct(pnode->Coord(),tmp2);
//   }
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tFinitElement::LocalBasis (const Tensor1s& lc_, Tensor1s& e1_,
//                                                          Tensor1s& e2_,
//                                                          Tensor1s& e3_) const
//{
// e1_.Assign0();   e2_.Assign0();   e3_.Assign0();
// Tensor1s nablaNi, r1, r2, r3;
// for (cardinal_t i=1; i<=HowManyCoordApproxNodes(); ++i)
//   {
//    CoordShapeGrad(i,lc_,nablaNi);
//    r1 = r2 = r3 = CoordApproxNode(i).Coord();
//    e1_ += (r1 *= nablaNi(1));
//    e2_ += (r2 *= nablaNi(2));
//    e3_ += (r3 *= nablaNi(3));
//   }
// return e1_;
//}
////------------------------------------------------------------------------------
//const Tensor2s& tFinitElement::TensorOfRotationToLocal (const Tensor1s& lc_) const
//{
// tMapRotToLoc::iterator p = RotToLoc_DB.find(make_pair(lc_(1),make_pair(lc_(2),lc_(3))));// modify to find with given precision!!!
// if (p == RotToLoc_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<tMapRotToLoc::iterator,bool> ins = RotToLoc_DB.insert(tMapRotToLoc::value_type(make_pair(lc_(1),make_pair(lc_(2),lc_(3))),Tensor2s()));
//    Assert(ins.second, "failed insertion in tFinitElement::TensorOfRotationToLocal");
//    p = ins.first;
//#else
//    p = RotToLoc_DB.insert(tMapRotToLoc::value_type(make_pair(lc_(1),make_pair(lc_(2),lc_(3))),Tensor2s())).first;
//#endif
//    TensorOfRotationToLocal(lc_,p->second);
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
///*Tensor2s& tFinitElement::DisplacementGradient (const Tensor1s& lc_, Tensor2s& result_) const
//{
// result_.Assign0();   Tensor1s tmp1;   Tensor2s tmp2;   const tNode* pnode;
// cardinal_t j=0;
// for (cardinal_t i=1; i<=HowManyDisplApproxNodes(); ++i)
//   {
//    pnode = &Node(i);
//    if (!pnode->DisplacementIs0())
//      result_ += (DisplShapeGrad(i,lc_,tmp1).DiadeProduct(pnode->Displacement(),tmp2));
//   }
// if (!result_.is0())
//   result_.ScalarMultiply_left(JacobyMatrix_l2g(lc_,tmp2).Invert());
// return result_;
//}*/
////------------------------------------------------------------------------------
////real_t tFinitElement::Mass() const
////{
//// return /*Volume() */ Material().Density();
////}
////------------------------------------------------------------------------------
//Tensor2s& tFinitElement::NablaDispl_g(const Tensor1s& lc_, Tensor2s& result_) const
//{
// Tensor1s tmp1;  Tensor2s tmp2;
// bool not_added(true);
// for (cardinal_t i=1; i <= HowManyDisplApproxNodes(); ++i)
//   if (!DisplApproxNode(i).DisplacementIs0())
//    {
//      DisplShapeGrad(DisplApproxNode(i),lc_,tmp1);
////      tmp1.DotMultiply_left(JacobyMatrix_l2g(lc_,tmp2).Invert()); -- singular matrix because JacobyMatrix_l2g(lc_,tmp2) is non-virtual
//      tmp1.DotMultiply_left((tmp2 = JacobyMatrix_l2g(lc_)).Invert());
//      tmp1.DirectProduct(DisplApproxNode(i).Displacement(),tmp2);
//      if (not_added)
//            {
//             result_ = tmp2;
//             not_added = false;
//            }
//       else
//            result_ += tmp2;
//   }
// if (not_added)
//       return  result_.Assign0();
//
// return result_;
//}
////------------------------------------------------------------------------------
////--------------------------tNonIsoparametricFE---------------------------------
////------------------------------------------------------------------------------
//tFinitElement& tNonIsoparametricFE::DefineNextNode (const tNode& node_, bool coordApproxing_, bool displApproxing_)
//{
//#ifdef STRONGCHECK
// Assert(coordApproxing_ || displApproxing_, "Node approxes neither coord nor displ");
//#endif
// tFinitElement::DefineNextNode (node_);
// const tNode* null = NULL;   vector<tNodeAnd2FunsPtrs>::iterator ppapproxdata;
// if (coordApproxing_)
//   {
//    ppapproxdata = find_if(CoordApproxData.begin(),CoordApproxData.end(),tNodeAnd2FunsPtrs::Node_address_is(null));
//#ifdef STRONGCHECK
//    Assert(ppapproxdata<CoordApproxData.end(), "Attempt to link redundant coord approxing node to FE");
//#endif
//    ppapproxdata->pNode = &node_;
//   }
// if (displApproxing_)
//   {
//    ppapproxdata = find_if(DisplApproxData.begin(),DisplApproxData.end(),tNodeAnd2FunsPtrs::Node_address_is(null));
//#ifdef STRONGCHECK
//    Assert(ppapproxdata<DisplApproxData.end(), "Attempt to link redundant coord approxing node to FE");
//#endif
//    ppapproxdata->pNode = &node_;
//   }
// return *this;
//}
////------------------------------------------------------------------------------
//const tNode& tNonIsoparametricFE::CoordApproxNode (cardinal_t nodeNo_) const
//{
//#ifndef STRONGCHECK
// Assert(nodeNo_> 0 &&  nodeNo_<= CoordApproxData.size(), "Invalid coord approx node # in FE");
//#endif
// return *(CoordApproxData[--nodeNo_].pNode);
//}
////------------------------------------------------------------------------------
//const tNode& tNonIsoparametricFE::DisplApproxNode (cardinal_t nodeNo_) const
//{
//#ifndef STRONGCHECK
// Assert(nodeNo_> 0  &&  nodeNo_<= DisplApproxData.size(), "Invalid displ approx node # in FE");
//#endif
// return *(DisplApproxData[--nodeNo_].pNode);
//}
////------------------------------------------------------------------------------
////--------------------------tIsoparametricFE------------------------------------
////------------------------------------------------------------------------------
//template <typename T>
//T tIsoparametricFE::tArray<T>::operator()(cardinal_t i_) const
//{
//#ifdef STRONGCHECK
// Assert(i_>0 && i_<=this->size(), "invalid index of shape function or its grad in operator()");
//#endif
//
// return this->operator[](--i_);
//}
////------------------------------------------------------------------------------
//template <typename T>
//T tIsoparametricFE::tArray<T>::operator()(const tFinitElement* pFE_, const tNode* pnode_) const
//{
//#ifdef STRONGCHECK
// Assert(pFE_->HowManyNodes() == this->size(), "number of nodes in FE is not equal to number of shape funs(grads)");
//#endif
//
// typename vector<T>::const_iterator ptr = this->begin();
// vector<const tNode*>::const_iterator ppnode = pFE_->p1stNode();
// while (ptr != this->end()  &&  *ppnode != pnode_)
// {++ptr; ++ppnode;}
//#ifdef STRONGCHECK
// Assert(*ppnode == pnode_, "the node is not belongs to the FE");
//#endif
// return *ptr;
//}
////------------------------------------------------------------------------------
//#ifdef STRONGCHECK
//tFinitElement& tIsoparametricFE::DefineNextNode (const tNode& node_, bool coordApproxing_, bool displApproxing_)
//{
// Assert(coordApproxing_ && displApproxing_, "Node  does not approx coord or displ in isoparametric FE");
// tFinitElement::DefineNextNode (node_);
// return *this;
//}
//#endif
////------------------------------------------------------------------------------
////-----------------------------t1D_FE-------------------------------------------
////------------------------------------------------------------------------------
//const fIntegrate<1,GAUSS_ORDER_1D> t1D_FE::Integrator;
////------------------------------------------------------------------------------
//Tensor1s& t1D_FE::GlobalCoord (const Tensor1s& lc_, Tensor1s& result_) const
//{
// tFinitElement::GlobalCoord(lc_,result_);
//// if (!lc_.is0(2)  ||  !lc_.is0(3))
////   {
//    Tensor1s normal1, normal2, tmp;   LocalBasis(lc_,tmp,normal1,normal2);
//    (result_ += (normal1 *= lc_(2))) += (normal2 *= lc_(3));
////   }
// return result_;
//}
////------------------------------------------------------------------------------
//const Tensor2s& t1D_FE::JacobyMatrix_l2g (const Tensor1s& lc_) const
//{
// map<real_t,Tensor2s>::iterator p = JacobyMatrix_DB.find(lc_(1));
// if (p == JacobyMatrix_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<map<real_t,Tensor2s>::iterator,bool> ins = JacobyMatrix_DB.insert(map<real_t,Tensor2s>::value_type(lc_(1),Tensor2s()));
//    Assert(ins.second, "failed insertion of ElasTensor in t1D_FE::JacobyMatrix_l2g(const Tensor1s&)");
//    p = ins.first;
//#else
//    p = JacobyMatrix_DB.insert(map<real_t,Tensor2s>::value_type(lc_(1),Tensor2s())).first;
//#endif
//    tFinitElement::JacobyMatrix_l2g(lc_,p->second);
//    Tensor1s normal1, normal2, tmp;
//    LocalBasis(lc_,tmp,normal1,normal2);
//    // p->second = Second_vector_of_basis.DiadeProdact(normal1):
//    p->second(2,1) = normal1(1);
//    p->second(2,2) = normal1(2);
//    p->second(2,3) = normal1(3);
//    // p->second = Third_vector_of_basis.DiadeProdact(normal2):
//    p->second(3,1) = normal2(1);
//    p->second(3,2) = normal2(2);
//    p->second(3,3) = normal2(3);
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
//Tensor1s& t1D_FE::Tangent (const Tensor1s& lc_, Tensor1s& result_) const
//{
// result_.Assign0();
// Tensor1s nodalRadiusVector, shapeGrad;
// for (cardinal_t i=1; i<=HowManyCoordApproxNodes(); ++i)
//   {
//    CoordShapeGrad(i,lc_,shapeGrad);
//    nodalRadiusVector = Node(i).Coord();
//    result_+= (nodalRadiusVector *= shapeGrad(1));
//   }
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& t1D_FE::LocalBasis (const Tensor1s& lc_, Tensor1s& tangent_,
//                                                   Tensor1s& normal1_,
//                                                   Tensor1s& normal2_) const
//{
// Tangent(lc_,tangent_);
// cardinal_t i = fabs(tangent_(1))<fabs(tangent_(2))? 2:1;
// if (fabs(tangent_(i))<fabs(tangent_(3))) i=1;
//    else                                  ++i;
// tangent_.CrossProductWithBasisVector(i,normal2_).Normalize();
// normal2_.CrossProduct(tangent_,normal1_).Normalize();
// return tangent_;
//}
////------------------------------------------------------------------------------
//Tensor2s& t1D_FE::TensorOfRotationToLocal (const Tensor1s& lc_, Tensor2s& result_) const
//{
// Tensor1s rotatedBasis1, rotatedBasis2, rotatedBasis3;
// //Assigning: rB1 - tangent to the FE direction, rB2,rB3 - unit normals to the direction:
// LocalBasis(lc_,rotatedBasis1,rotatedBasis2,rotatedBasis3).Normalize();//tangent is normalised
// //now rB123 correspond to the basis of global coord system, rotated to the FE surface
// return result_.AssignRows(rotatedBasis1,rotatedBasis2,rotatedBasis3);
//}
////------------------------------------------------------------------------------
//const SymmetricTensor4s& t1D_FE::MaterialElasTensor (const Tensor1s& lc_) const
//{
// if (StateIs1dStress)
//      return pMaterial->ElasTensor_1D(tFinitElement::TensorOfRotationToLocal(lc_));
//  else
//      return tFinitElement::MaterialElasTensor(lc_);
//}
////------------------------------------------------------------------------------
//real_t t1D_FE::Jacobian_l2g (const Tensor1s& lc_) const
//{
// map<real_t,real_t>::iterator p = Jacobian_DB.find(lc_(1));
// if (p == Jacobian_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<map<real_t,real_t>::iterator,bool> ins = Jacobian_DB.insert(map<real_t,real_t>::value_type(lc_(1),JacobyMatrix_l2g(lc_).Det()));
//    Assert(ins.second, "failed insertion of ElasTensor in t1D_FE::Jacobian_l2g(const Tensor1s&)");
//    p = ins.first;
//#else
//    p = Jacobian_DB.insert(map<real_t,real_t>::value_type(lc_(1),JacobyMatrix_l2g(lc_).Det())).first;
//#endif
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
//real_t t2D_FE::Jacobian_l2g (const Tensor1s& lc_) const
//{
// map<pair<real_t,real_t>,real_t>::iterator p = Jacobian_DB.find(make_pair(lc_(1),lc_(2)));
// if (p == Jacobian_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<map<pair<real_t,real_t>,real_t>::iterator,bool> ins = Jacobian_DB.insert(map<pair<real_t,real_t>,real_t>::value_type(make_pair(lc_(1),lc_(2)),JacobyMatrix_l2g(lc_).Det()));
//    Assert(ins.second, "failed insertion of ElasTensor in t2D_FE::Jacobian_l2g(const Tensor1s&)");
//    p = ins.first;
//#else
//    p = Jacobian_DB.insert(map<pair<real_t,real_t>,real_t>::value_type(make_pair(lc_(1),lc_(2)),JacobyMatrix_l2g(lc_).Det())).first;
//#endif
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
//real_t t3D_FE::Jacobian_l2g (const Tensor1s& lc_) const
//{
// map<pair<pair<real_t,real_t>,real_t>,real_t>::iterator p = Jacobian_DB.find(make_pair(make_pair(lc_(1),lc_(2)),lc_(3)));
// if (p == Jacobian_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<map<pair<pair<real_t,real_t>,real_t>,real_t>::iterator,bool> ins
//    = Jacobian_DB.insert(map<pair<pair<real_t,real_t>,real_t>,real_t>::value_type(make_pair(make_pair(lc_(1),lc_(2)),lc_(3)),JacobyMatrix_l2g(lc_).Det()));
//    Assert(ins.second, "failed insertion of ElasTensor in t3D_FE::Jacobian_l2g(const Tensor1s&)");
//    p = ins.first;
//#else
//    p = Jacobian_DB.insert(map<pair<pair<real_t,real_t>,real_t>,real_t>::value_type(make_pair(make_pair(lc_(1),lc_(2)),lc_(3)),JacobyMatrix_l2g(lc_).Det())).first;
//#endif
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
////-----------------------------t2D_FE-------------------------------------------
////------------------------------------------------------------------------------
//Tensor1s& t2D_FE::LocalBasis (const Tensor1s& lc_, Tensor1s& tangent1_,
//                                                          Tensor1s& tangent2_,
//                                                          Tensor1s& normal_) const
//{
// tangent1_.Assign0();   tangent2_.Assign0();
// Tensor1s nablaNi, radiusVectOfNodei;
//
// for (cardinal_t i=1; i<=HowManyCoordApproxNodes(); ++i)
//   {
//    CoordShapeGrad(i,lc_,nablaNi);
//    radiusVectOfNodei = normal_ = CoordApproxNode(i).Coord();//here normal_ is used as tmp storage
//    tangent1_ += (radiusVectOfNodei *= nablaNi(1));
//    tangent2_ += (normal_           *= nablaNi(2));
//   }
// tangent1_.CrossProduct(tangent2_,normal_).Normalize();
// return tangent1_;
//}
////------------------------------------------------------------------------------
//Tensor1s& t2D_FE::GlobalCoord (const Tensor1s& lc_, Tensor1s& result_) const
//{
// tFinitElement::GlobalCoord(lc_,result_);
// Tensor1s tmp;
//// if (!lc_.is0(3))
//                   result_ += (UnitNormal(lc_,tmp) *= lc_(3));
// return result_;
//}
////------------------------------------------------------------------------------
//const Tensor2s& t2D_FE::JacobyMatrix_l2g (const Tensor1s& lc_) const
//{
// map<pair<real_t,real_t>,Tensor2s>::iterator p = JacobyMatrix_DB.find(make_pair(lc_(1),lc_(2)));
// if (p == JacobyMatrix_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<map<pair<real_t,real_t>,Tensor2s>::iterator,bool> ins = JacobyMatrix_DB.insert(map<pair<real_t,real_t>,Tensor2s>::value_type(make_pair(lc_(1),lc_(2)),Tensor2s()));
//    Assert(ins.second, "failed insertion of ElasTensor in t2D_FE::JacobyMatrix_l2g(const Tensor1s&)");
//    p = ins.first;
//#else
//    p = JacobyMatrix_DB.insert(map<pair<real_t,real_t>,Tensor2s>::value_type(make_pair(lc_(1),lc_(2)),Tensor2s())).first;
//#endif
//    tFinitElement::JacobyMatrix_l2g(lc_,p->second);
//    Tensor1s normal1;
//    UnitNormal(lc_,normal1);
//    /* p->second += ThirdBasicVector.DiadeProdact(normal1):*/
//    p->second(3,1) += normal1(1);
//    p->second(3,2) += normal1(2);
//    p->second(3,3) += normal1(3);
//    // if (!lc_.is0(3)) {/*must be defined for use of this function beyond the median surface of FE*/}
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
//Tensor2s& t2D_FE::TensorOfRotationToLocal (const Tensor1s& lc_, Tensor2s& result_) const
//{
// Tensor1s rotatedBasis1, rotatedBasis2, rotatedBasis3;
// //Assigning: rB1, rB2 - coord.vectors of FE surface, rB3 - normal to the surface:
// LocalBasis(lc_,rotatedBasis1,rotatedBasis2,rotatedBasis3);
// //Because rB1, rB2 are not orthogonal, rB2 is replaced by vector orthogonal to rB1,rB3:
// rotatedBasis1.Normalize();
// rotatedBasis3.CrossProduct(rotatedBasis1,rotatedBasis2);
// //now rB123 corresponds to the basis of global coord system, rotated to FE surface
// return result_.AssignRows(rotatedBasis1,rotatedBasis2,rotatedBasis3);
//}
////------------------------------------------------------------------------------
//const SymmetricTensor4s& t2D_FE::MaterialElasTensor (const Tensor1s& lc_) const
//{
// if (StateIs2dStress)
//      return pMaterial->ElasTensor_2D(tFinitElement::TensorOfRotationToLocal(lc_));
//  else
//      return tFinitElement::MaterialElasTensor(lc_);
//}
////------------------------------------------------------------------------------
////-----------------------------t3D_FE-------------------------------------------
////------------------------------------------------------------------------------
//const Tensor2s& t3D_FE::JacobyMatrix_l2g (const Tensor1s& lc_) const
//{
// map<pair<pair<real_t,real_t>,real_t>,Tensor2s>::iterator p = JacobyMatrix_DB.find(make_pair(make_pair(lc_(1),lc_(2)),lc_(3)));
// if (p == JacobyMatrix_DB.end())
//   {
//#ifdef STRONGCHECK
//    pair<map<pair<pair<real_t,real_t>,real_t>,Tensor2s>::iterator,bool> ins
//     = JacobyMatrix_DB.insert(map<pair<pair<real_t,real_t>,real_t>,Tensor2s>::value_type(make_pair(make_pair(lc_(1),lc_(2)),lc_(3)),Tensor2s()));
//    Assert(ins.second, "failed insertion of ElasTensor in t3D_FE::JacobyMatrix_l2g(const Tensor1s&)");
//    p = ins.first;
//#else
//    p = JacobyMatrix_DB.insert(map<pair<pair<real_t,real_t>,real_t>,Tensor2s>::value_type(make_pair(make_pair(lc_(1),lc_(2)),lc_(3)),Tensor2s())).first;
//#endif
//    tFinitElement::JacobyMatrix_l2g(lc_,p->second);
//   }
// return p->second;
//}
////------------------------------------------------------------------------------
//Tensor2s& t3D_FE::TensorOfRotationToLocal (const Tensor1s& lc_, Tensor2s& result_) const
//{
// Tensor1s rotatedBasis1, rotatedBasis2, rotatedBasis3;
// LocalBasis(lc_,rotatedBasis1,rotatedBasis2,rotatedBasis3);
// rotatedBasis1.Normalize().CrossProduct(rotatedBasis2,rotatedBasis3);
// rotatedBasis3.Normalize().CrossProduct(rotatedBasis1,rotatedBasis2);
// return result_.AssignRows(rotatedBasis1,rotatedBasis2,rotatedBasis3);
//}
////------------------------------------------------------------------------------
////-----------------------------tRectangle---------------------------------------
////------------------------------------------------------------------------------
//const fIntegrate<2,GAUSS_ORDER_2D> tRectangle::Integrator;
////------------------------------------------------------------------------------
////-----------------------------tIsoRod2ConstSec--------------------------
////------------------------------------------------------------------------------
//real_t tIsoRod2ConstSec::ShapeFun1(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// return 0.5l*(1.l - lc_(1));
//}
////------------------------------------------------------------------------------
//real_t tIsoRod2ConstSec::ShapeFun2(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// return 0.5l*(1.l + lc_(1));
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoRod2ConstSec::ShapeGrad1(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1,-0.5l);
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoRod2ConstSec::ShapeGrad2(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1,0.5l);
// return result_;
//}
////------------------------------------------------------------------------------
//const tIsoRod2ConstSec::tPtrShapeFunArray     tIsoRod2ConstSec::pShapeFunctions;
//const tIsoRod2ConstSec::tPtrShapeFunGradArray tIsoRod2ConstSec::pShapeFunGrads;
////------------------------------------------------------------------------------
//real_t tIsoRod2ConstSec::CoordShapeFun (cardinal_t no_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(no_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoRod2ConstSec::CoordShapeGrad (cardinal_t no_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(no_)))(lc_,result_);
//}
////------------------------------------------------------------------------------
///*real_t   tIsoRod2ConstSec::DisplShapeFun (cardinal_t no_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(no_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoRod2ConstSec::DisplShapeGrad(cardinal_t no_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(no_)))(lc_,result_);
//}*/
////------------------------------------------------------------------------------
//real_t tIsoRod2ConstSec::CoordShapeFun (const tNode& node_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(this,&node_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoRod2ConstSec::CoordShapeGrad (const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
//}
////------------------------------------------------------------------------------
///*real_t tIsoRod2ConstSec::DisplShapeFun (const tNode& node_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(this,&node_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1& tIsoRod2ConstSec::DisplShapeGrad (const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
//}*/
////------------------------------------------------------------------------------
////-----------------------------tIsoQuad4ConsThick-------------------------------
////------------------------------------------------------------------------------
//real_t tIsoQuad4ConsThick::ShapeFun1(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// return 0.25l*(1.l + lc_(1))*(1.l - lc_(2));
//}
////------------------------------------------------------------------------------
//real_t tIsoQuad4ConsThick::ShapeFun2(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// return 0.25l*(1.l - lc_(1))*(1.l - lc_(2));
//}
////------------------------------------------------------------------------------
//real_t tIsoQuad4ConsThick::ShapeFun3(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// return 0.25l*(1.l - lc_(1))*(1.l + lc_(2));
//}
////------------------------------------------------------------------------------
//real_t tIsoQuad4ConsThick::ShapeFun4(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// return 0.25l*(1.l + lc_(1))*(1.l + lc_(2));
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::ShapeGrad1(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();//to do assign only last component
// result_.Assign(1, 0.25l*( 1.l - lc_(2)));
// result_.Assign(2, 0.25l*(-1.l - lc_(1)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::ShapeGrad2(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, 0.25l*(-1.l + lc_(2)));
// result_.Assign(2, 0.25l*(-1.l + lc_(1)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::ShapeGrad3(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, tMessage("value of local coord must be between -1 & 1\nvaluse are:\n")<<lc_(1)<<'\n'<<lc_(2));
//#endif
// result_.Assign0();
// result_.Assign(1, 0.25l*(-1.l - lc_(2)));
// result_.Assign(2, 0.25l*( 1.l - lc_(1)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::ShapeGrad4(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l, "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, 0.25l*( 1.l + lc_(2)));
// result_.Assign(2, 0.25l*( 1.l + lc_(1)));
// return result_;
//}
////------------------------------------------------------------------------------
//tIsoQuad4ConsThick::tPtrShapeFunArray::tPtrShapeFunArray():
//    tIsoparametricFE::tArrayOfPtrsToShapeFunctions(4)
//{
// operator[](0) = ShapeFun1;
// operator[](1) = ShapeFun2;
// operator[](2) = ShapeFun3;
// operator[](3) = ShapeFun4;
//}
////------------------------------------------------------------------------------
//const tIsoQuad4ConsThick::tPtrShapeFunArray tIsoQuad4ConsThick::pShapeFunctions;
////------------------------------------------------------------------------------
//tIsoQuad4ConsThick::tPtrShapeFunGradArray::tPtrShapeFunGradArray():
//    tIsoparametricFE::tArrayOfPtrsToShapeFunGradients(4)
//{
// operator[](0) = ShapeGrad1;
// operator[](1) = ShapeGrad2;
// operator[](2) = ShapeGrad3;
// operator[](3) = ShapeGrad4;
//}
////------------------------------------------------------------------------------
//const tIsoQuad4ConsThick::tPtrShapeFunGradArray tIsoQuad4ConsThick::pShapeFunGrads;
////------------------------------------------------------------------------------
//real_t tIsoQuad4ConsThick::CoordShapeFun (cardinal_t no_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(no_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::CoordShapeGrad (cardinal_t no_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(no_)))(lc_,result_);
//}
////------------------------------------------------------------------------------
///*real_t   tIsoQuad4ConsThick::DisplShapeFun (cardinal_t no_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(no_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::DisplShapeGrad(cardinal_t no_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(no_)))(lc_,result_);
//}*/
////------------------------------------------------------------------------------
//real_t tIsoQuad4ConsThick::CoordShapeFun (const tNode& node_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(this,&node_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::CoordShapeGrad (const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
//}
////------------------------------------------------------------------------------
///*real_t tIsoQuad4ConsThick::DisplShapeFun (const tNode& node_, const Tensor1& lc_) const
//{
// return (*(pShapeFunctions(this,&node_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoQuad4ConsThick::DisplShapeGrad (const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
//}*/
////------------------------------------------------------------------------------
////------------------------------------------------------------------------------
////------------------------------------------------------------------------------
////------------------------------------------------------------------------------
//const fIntegrate<3,GAUSS_ORDER_3D> tIsoParallelepiped8::Integrator;
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun1(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l+lc_(1))*(1.l-lc_(2))*(1.l-lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun2(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l-lc_(1))*(1.l-lc_(2))*(1.l-lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun3(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//        "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l-lc_(3))*(1.l+lc_(2))*(1.l-lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun4(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l+lc_(1))*(1.l+lc_(2))*(1.l-lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun5(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l+lc_(1))*(1.l-lc_(2))*(1.l+lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun6(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l-lc_(1))*(1.l-lc_(2))*(1.l+lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun7(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l-lc_(1))*(1.l+lc_(2))*(1.l+lc_(3));
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::ShapeFun8(const Tensor1s& lc_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// return 0.125l*(1.l+lc_(1))*(1.l+lc_(2))*(1.l+lc_(3));
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad1(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1,  0.125l*(1.l-lc_(2))*(1.l-lc_(3)));
// result_.Assign(2, -0.125l*(1.l+lc_(1))*(1.l-lc_(3)));
// result_.Assign(3, -0.125l*(1.l+lc_(1))*(1.l-lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad2(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, -0.125l*(1.l-lc_(2))*(1.l-lc_(3)));
// result_.Assign(2, -0.125l*(1.l-lc_(1))*(1.l-lc_(3)));
// result_.Assign(3, -0.125l*(1.l-lc_(1))*(1.l-lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad3(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, -0.125l*(1.l+lc_(2))*(1.l-lc_(3)));
// result_.Assign(2,  0.125l*(1.l-lc_(1))*(1.l-lc_(3)));
// result_.Assign(3, -0.125l*(1.l-lc_(1))*(1.l+lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad4(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1,  0.125l*(1.l+lc_(2))*(1.l-lc_(3)));
// result_.Assign(2,  0.125l*(1.l+lc_(1))*(1.l-lc_(3)));
// result_.Assign(3, -0.125l*(1.l+lc_(1))*(1.l+lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad5(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1,  0.125l*(1.l-lc_(2))*(1.l+lc_(3)));
// result_.Assign(2, -0.125l*(1.l+lc_(1))*(1.l+lc_(3)));
// result_.Assign(3,  0.125l*(1.l+lc_(1))*(1.l-lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad6(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, -0.125l*(1.l-lc_(2))*(1.l+lc_(3)));
// result_.Assign(2, -0.125l*(1.l-lc_(1))*(1.l+lc_(3)));
// result_.Assign(3,  0.125l*(1.l-lc_(1))*(1.l-lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad7(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, -0.125l*(1.l+lc_(2))*(1.l+lc_(3)));
// result_.Assign(2,  0.125l*(1.l-lc_(1))*(1.l+lc_(3)));
// result_.Assign(3,  0.125l*(1.l-lc_(1))*(1.l+lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::ShapeGrad8(const Tensor1s& lc_, Tensor1s& result_)
//{
//#ifdef STRONGCHECK
// Assert(lc_(1)>=-1.l && lc_(1)<=1.l && lc_(2)>=-1.l && lc_(2)<=1.l && lc_(3)>=-1.l && lc_(3)<=1.l,
//  "value of local coord must be between -1 & 1");
//#endif
// result_.Assign0();
// result_.Assign(1, 0.125l*(1.l+lc_(2))*(1.l+lc_(3)));
// result_.Assign(2, 0.125l*(1.l+lc_(1))*(1.l+lc_(3)));
// result_.Assign(3, 0.125l*(1.l+lc_(1))*(1.l+lc_(2)));
// return result_;
//}
////------------------------------------------------------------------------------
//tIsoParallelepiped8::tPtrShapeFunArray::tPtrShapeFunArray():
//        tIsoparametricFE::tArrayOfPtrsToShapeFunctions(8)
//{
//  operator[](0) = ShapeFun1;
//  operator[](1) = ShapeFun2;
//  operator[](2) = ShapeFun3;
//  operator[](3) = ShapeFun4;
//  operator[](4) = ShapeFun5;
//  operator[](5) = ShapeFun6;
//  operator[](6) = ShapeFun7;
//  operator[](7) = ShapeFun8;
//}
////------------------------------------------------------------------------------
//const tIsoParallelepiped8::tPtrShapeFunArray tIsoParallelepiped8::pShapeFunctions;
////------------------------------------------------------------------------------
//tIsoParallelepiped8::tPtrShapeFunGradArray::tPtrShapeFunGradArray():
//        tIsoparametricFE::tArrayOfPtrsToShapeFunGradients(8)
//{
//  operator[](0) = ShapeGrad1;
//  operator[](1) = ShapeGrad2;
//  operator[](2) = ShapeGrad3;
//  operator[](3) = ShapeGrad4;
//  operator[](4) = ShapeGrad5;
//  operator[](5) = ShapeGrad6;
//  operator[](6) = ShapeGrad7;
//  operator[](7) = ShapeGrad8;
//}
////------------------------------------------------------------------------------
//const tIsoParallelepiped8::tPtrShapeFunGradArray tIsoParallelepiped8::pShapeFunGrads;
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::CoordShapeFun (cardinal_t no_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(no_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::CoordShapeGrad (cardinal_t no_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(no_)))(lc_,result_);
//}
////------------------------------------------------------------------------------
//real_t tIsoParallelepiped8::CoordShapeFun (const tNode& node_, const Tensor1s& lc_) const
//{
// return (*(pShapeFunctions(this,&node_)))(lc_);
//}
////------------------------------------------------------------------------------
//Tensor1s& tIsoParallelepiped8::CoordShapeGrad (const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
//{
// return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
//}
////------------------------------------------------------------------------------
////template<typename T>
////T& t3D_FE::Integrate3D (const tFinitElement::fIntegrand<T>& function_, T& result_) const
////{
////  Tensor1s localCoord(0.l);
////  real_t gaussWeightX, gaussWeightY, gaussWeightZ;   cardinal_t i=0, j, k;
////  const vector<tFinitElement::tGaussCoeff> &gausscoeffsX = *pGaussCoeffsX,
////                                           &gausscoeffsY = *pGaussCoeffsY,
////                                           &gausscoeffsZ = *pGaussCoeffsZ;
////  for (; i<gausscoeffsX.size(); ++i)
////    {
////     gausscoeffsX[i].WhatCoeffs(localCoord(1),gaussWeightX);
////     for (j=0; j<gausscoeffsY.size(); ++j)
////       {
////        gausscoeffsY[j].WhatCoeffs(localCoord(2),gaussWeightY);
////        for (k=0; k<gausscoeffsZ.size(); ++k)
////          {
////           gausscoeffsZ[k].WhatCoeffs(localCoord(3),gaussWeightZ);
////           result_ += (function_(localCoord) *=
////                       (gaussWeightZ * gaussWeightY * gaussWeightX *
////                        Jacobian_l2g(localCoord)));
////          }
////       }
////    }
////  return result_;
////}
