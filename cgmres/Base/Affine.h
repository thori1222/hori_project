/*
 *  Copyright (c) 2003-2008, Shoichi Hasegawa and Springhead development team 
 *  All rights reserved.
 *  This software is free software. You can freely use, distribute and modify this 
 *  software. Please deal with this software under one of the following licenses: 
 *  This license itself, Boost Software License, The MIT License, The BSD License.   
 */
#ifndef AFFINE_H
#define AFFINE_H

/** \addtogroup gpLinearAlgebra */
//@{
/** \defgroup gpAffine @!JAAffine�s��N���X@!ENA Affine Matrix@!*
	\section secAffineIntro @!JA�͂��߂�@!ENIntroduction@!*
		@!JA
		���̃y�[�W��Affine�s��N���X(Spr::TAffine, Spr::TAffine2)�̐����ł��D
		Affine�s��N���X���C�u�����́C3D�V�~�����[�V�����ɕK�{��
		Affine�s���C++�̃N���X�ɂ������̂ł��D
		@!EN
		This is a document for Affine Matrix Library.
		Affine Matrix Library is a set of classes for affine matrix,
		which are necessary for 3D simulation.
		@!*
	\section secAffineUsage @!JA�g����@!ENUsage@!*
	@!JA
	Affine �s��N���X���C�u�����́C�w�b�_�t�@�C����������Ȃ�
	�N���X���C�u�����Ȃ̂�, TAffine.h, TinyVec.h, TinyMat.h, TMatrix.h, TMatrixUtility.h, TVector.h
	�𓯂��t�H���_�ɓ���Ă����C.cpp�t�@�C������w�b�_���C���N���[�h���邾����
	�g�p�ł��܂��D
	@!EN
	Affine Matrix Library consist only from some header files.
	You can use this library just coping some header files
	(TAffine.h, TinyVec.h, TinyMat.h, TMatrix.h, TMatrixUtility.h, TVector.h)
	and including them from your .cpp files.
	@!*
	\subsection secAffineSample @!JA �T���v�� @!EN sample@!*
	\verbatim
#include "Affine.h"                         //  TAffine�s�񃉃C�u�����̃C���N���[�h����D
#include <iostream>

using namespace Spr;	//  @!JA Affine�s��N���X��Spr���O��Ԃ̒��Ő錾����Ă���D@!EN Affine matrix class is declared in the name space of Spr.@!*

void main(){
    Affinef af=Affinef::Rad(Rad(30), 'z');  //  @!JA�v�f��float��TAffine�s���錾. 
                                                @!ENDeclare an Affine matrix of float element.
                                                @!*
                                                @!JAz�����30�x��]�s��ɏ�����
                                                @!ENInitializing as a matrix of rotation of 30 degrees around z axis.
                                                @!*
                                                
    Vec3f vec(1,0,0);                       //  @!JA�v�f��float��3�����̃x�N�g����錾
                                                @!ENDeclare a vector of float element.
                                                @!*
    std::cout << af;
    std::cout << vec << std::endl;
    std::cout << af * vec << std::endl;
}\endverbatim
	\subsection secAffineVecFunc @!JA Affine�s��E�x�N�g���̃����o�Ɖ��Z  @!EN Functions of Affine matrixs and vectors. @!*
	@!JA
	���ʂɉ��Z���ł��܂��D
	<ul>
	<li> +:�a, -:��, *:��/����, /:�萔����1
	<li> ==:��r, =:���
	<li> <<:�o��, >>:����
	<li> %:�x�N�g���̊O��
	</ul>
	@!EN Normal calculations are supported.
	<ul>
	<li> +:addition, -:subtraction, *:multiplication/inner product, /:division
	<li> ==:comparison, =:assign
	<li> <<:output, >>:input
	<li> %:cross product
	</ul>
	@!*
	@!JA
	Affine�ϊ��́C
	\verbatim
	TAffine<float> af; TVec3<float> v, af_v;
	af_v = af * v;\endverbatim
	�Ƃ���΂ł��܂��D
	@!EN
	Affine transformation is executed by:
	\verbatim
	TAffine<float> af; TVec3<float> v, af_v;
	af_v = af * v;\endverbatim
	@!*
	@!JA
	�܂��C���̂悤��TAffine�s��̕��������o�����Ƃ��ł��܂��D
	<ul>
		<li> af.Ex():	X�����x�N�g���D(3�����x�N�g��)�D
		<li> af.Ey():	Y�����x�N�g���D(3�����x�N�g��)�D
		<li> af.Ez():	Z�����x�N�g���D(3�����x�N�g��)�D
		<li> af.Trn():	���s�ړ������D(3�����x�N�g��)�D
		<li> af.Rot():	��]�ϊ������D(3�~�R�s��)�D
	</ul>
	@!EN
	You can extract some portion of the Affine matrix as:
	<ul>
		<li> af.Ex():	base vector for x axis (3 dimensional vector)
		<li> af.Ey():	base vector for y axis (3 dimensional vector)
		<li> af.Ez():	base vector for z axis (3 dimensional vector)
		<li> af.Trn():	portion of translation�D(3 dimensional vector)
		<li> af.Rot():	portion of rotation�D(3x3 matrix)
	</ul>
	@!*
	@!JA �����ւ̑���Ȃǂ��ł��܂��D
	@!EN You can assign portion of affine matrx.
	@!*
	\verbatim
	TAffine<float> af;
	af.Pos() = Vec3f(10,0,0);
	af.Rot() = TAffine<float>::Rot(Rad(30), 'x').Rot() * af.Rot();\endverbatim
	@!JA �x�N�g���͎��̂悤�ȃ����o�֐��������܂��D
	@!EN a vector has following members.
	@!*
	<ul>
	<li> unit(): @!JA �������������P�ʃx�N�g����Ԃ��D
				 @!EN return a unit vector.
				 @!*
	<li> norm(): @!JA �x�N�g���̑傫��(�m����)��Ԃ��D
				 @!EN return the size (norm) of a vector.
				 @!*
	</ul>
	\subsection secAffineConstruct Affine�s��̏����� 
	TAffine�s��(Spr::TAffine)�ɂ͕֗��ȃR���X�g���N�^�⏉�����֐���p�ӂ��܂����D
	�������������ŏЉ�܂��D
	<ul>
		<li> Spr::TAffine::Trn (T px, T py, T pz):	
				���s�ړ�����s��ɏ������DT��TAffine<T>��T�Dfloat��double�Ȃǂŗǂ��D
		<li> Spr::TAffine::Rot(element_type th, char axis):
				��]�s���Ԃ��Dth�̓��W�A���Daxis�́C'x', 'y', 'z'�Delement_type��T�̂��ƁD
		<li> Spr::TAffine::ProjectionD3D (TVec3 screen, TVec2 size, T front=1.0f, T back=10000.0f):
				D3D�p�ˉe�s��Ƃ��ď������D
		<li> Spr::TAffine::ProjectionGL (TVec3 screen, TVec2 size, T front=1.0f, T back=10000.0f):
		<br>
			OpenGL�p�ˉe�s��Ƃ��ď�����(-Z���O)�D
			<ul>
				<li> screen  �J�������猩���r���[�|�[�g�̒��S�̈ʒu
				<li> size    �r���[�|�[�g�̃T�C�Y
				<li> front   ��O�̃N���b�s���O���ʂƃJ�����̋���
				<li> back    ���̃N���b�s���O���ʂƃJ�����̋���
			</ul>
		<li> Spr::TAffine::LookAtGL (TVec3 pos, TVec3 diry)
			�ʒu�͂��̂܂܂ŁCpos��-Ez(), diry �� Ey()�������悤��Affine�s��D
			OpenGL��gluLookAt�Ɠ����D
	</ul>
	
	\section thanks �ӎ�
	LU�����C�t�s��C�K�E�X�����@�Ȃǂ̍s��v�Z�A���S���Y���́C<br>
    �u�w�b����ɂ��ŐV�A���S���Y�����T�x�S�\�[�X�R�[�h�v<br>
    ftp://ftp.matsusaka-u.ac.jp/pub/algorithms<br>
	���� ���F Haruhiko Okumura<br>
	�����ς��ė��p�����Ă��������܂����D
	���R�ɃR�[�h���g����悤���J���Ă��������Ă��肪�Ƃ��������܂��D

	@section secRefLA ���t�@�����X
		@ref gpLinearAlgebra
*/
//@{


#include "TMatrix.h"
#include "TMatrixUtility.h"
#include "TinyVec.h"
#include "TinyMat.h"
#include "VectorDebug.h"

/**	@file Affine.h 2/3�����A�t�B���s��*/

#ifndef PTM_PACK	//	�P�̂Ŏg�p����ꍇ�́Cnamespace �ɓ���Ȃ�
namespace Spr {
#endif

#undef M_PI
#ifdef __BORLANDC__
#define M_PI 3.14159265358979323846
#else
///	�~������
const double M_PI = 3.14159265358979323846;
#endif

#undef abs
#ifdef __BORLANDC__
/*	��Βl�DBCB6���Ⴄ�֐��������N���Ă��܂�(���Ԃ�o�O)�̂ŁCtemplate �͎g�p���Ă��Ȃ�
	std::abs �ƂԂ��邽�߂��Ǝv����D	*/
#define DEF_ABS_FUNC(T)		inline T abs(T t){ return t > T()  ?  t  :  -t; }
DEF_ABS_FUNC(float)
DEF_ABS_FUNC(double)
DEF_ABS_FUNC(char)
DEF_ABS_FUNC(int)
DEF_ABS_FUNC(long)
#else
template <class T> T abs(T t){ return t > T()  ?  t  :  -t; }
#endif

#undef sign
///	����(���Ȃ�1, ���Ȃ� -1 ��Ԃ�)
template <class T> T sign(T t){
	return t > T()  ?  T(1)  :  T(-1);
}

#if _MSC_VER < 12

#undef min
/// ����������Ԃ�
template <class T> T min(T a, T b){return a < b ? a : b;}

#undef max
/// �傫������Ԃ�
template <class T> T max(T a, T b){return a > b ? a : b;}

#endif

/// ���ς��Ƃ�
template <class T> T ave(T a, T b){return T(0.5 * (a + b));}

///	�x�����W�A���ɕϊ�
inline double Rad(double deg){
	return ((double)deg/360*2*M_PI);
}
inline float Radf(double deg){
	return (float)((double)deg/360*2*M_PI);
}
///	���W�A����x�ɕϊ�
inline double Deg(double rad){
	return (rad/(2*M_PI)) * 360;
}
inline float Degf(double rad){
	return (float)(rad/(2*M_PI)) * 360;
}
///	2��
template <class SC>
inline SC Square(SC x){
	return x*x;
}
///	2x2�s��̍s��
template <class SC>
inline SC Det2(SC a, SC b, SC c, SC d){
	return ((a)*(d) - (b)*(c));
}

//-----------------------------------------------------------------------------
//	TAffine2

/**	TAffine2�s��(��],�g��,���s�ړ���\��)�s��.
	�T�v�́C\ref pgAffine �Q�ƁD
*/
template <class T>
class TAffine2:public PTM::TMatrixBase<3,3,
	PTM::TMatrixDescCol< TAffine2<T>, PTM::TMatrixRow<3,3,T>, 3,3,3,T> >{
public:
	typedef PTM::TMatrixDescCol< TAffine2<T>, PTM::TMatrixRow<3,3,T>, 3,3,3,T> desc;
	typedef PTM::TMatrixBase<3,3,desc> base_type;
	///	��{�I�ȃ����o�̒�` @see ::DEF_MATRIX_BASIC_MEMBER
	DEF_MATRIX_BASIC_MEMBER(TAffine2);
	union{
		struct{
			T xx, xy, xz;
			T yx, yy, yz;
			T px, py, pz;
		};
		T data[3][3];
	};
	///	�v�f�̃A�N�Z�X
	element_type& item_impl(size_t i, size_t j){ return data[j][i]; }
	const element_type& item_impl(size_t i, size_t j) const { return data[j][i]; }
	
	/**@name	���x�N�g���ւ̃A�N�Z�X	*/
	//@{
	/// 
	TVec2<element_type>& Ex(){
		return *(TVec2<element_type>*) &this->item(0,0);
	}
	/// 
	const TVec2<element_type>& Ex() const{
		return *(TVec2<element_type>*) &this->item(0,0);
	}
	/// 
	TVec2<element_type>& Ey(){
		return *(TVec2<element_type>*) &this->item(0,1);
	}
	/// 
	const TVec2<element_type>& Ey() const{
		return *(TVec2<element_type>*) &this->item(0,1);
	}
	/// 
	TVec2<element_type>& Trn(){
		return *(TVec2<element_type>*) &this->item(0,2);
	}
	/// 
	const TVec2<element_type>& Trn() const{
		return *(TVec2<element_type>*) &this->item(0,2);
	}
	/// 
	TVec2<element_type>& Pos(){ return Trn(); }
	/// 
	const TVec2<element_type>& Pos() const { return Trn(); }
	//@}

	/**@name	�v�f�ւ̃A�N�Z�X	*/
	//@{
	/// 
	element_type& ExX() {return Ex().X();}
	/// 
	element_type& ExY() {return Ex().Y();}
	/// 
	element_type& EyX() {return Ey().X();}
	/// 
	element_type& EyY() {return Ey().Y();}
	/// 
	element_type& TrnX() {return Trn().X();}
	/// 
	element_type& TrnY() {return Trn().Y();}
	///	TrnX()�̕ʖ�
	element_type& PosX() {return Trn().X();}
	///	TrnY()�̕ʖ�
	element_type& PosY() {return Trn().Y();}
	//@}

	///	��]�g��ϊ��������o��.
	PTM::TSubMatrixCol<2,2, desc>& Rot() { return sub_matrix(0,0,PTM::TSubMatrixCol<2,2, desc>()); }
	///	��]�g��ϊ��������o�� (const��).
	const PTM::TSubMatrixCol<2,2, desc>& Rot() const { return sub_matrix(0,0,PTM::TSubMatrixCol<2,2, desc>()); }
	
	/**@name	�������ƍ\�z	*/
	//@{
	///	�P�ʍs��
	static TAffine2<T> Unit(){
		TAffine2<T> y;
		PTM::init_unitize(y);
		return y;
	}
	///	���s�ړ�
	static TAffine2<T> Trn(element_type px, element_type py){
		TAffine2<T> y;
		y.Trn().X() = px;
		y.Trn().Y() = py;
		return y;
	}
	///	��]�C�Ȃ����������ЂƂ���VC.net�ŃG���[�ɂȂ�D
	static TAffine2<T> Rot(element_type th, int d=0){
		TAffine2 y;
		PTM::init_rot(y.Rot(), th);
		return y;
	}
	/// �g��
	static TAffine2<T> Scale(element_type sx, element_type sy){
		TAffine2<T> y;
		y.item(0, 0) = sx; y.item(1, 1) = sy;
		return y;
	}
	///�R���X�g���N�^
	void set_default(){PTM::init_unitize(*this);}
	//@}
};

///	TAffine2�ƃx�N�g���̊|���Z
template <class TD, class TV>
TVec2<TV> operator * (const PTM::TMatrixBase<3,3, TD>& a, const TVec2<TV>& b){
	TVec2<TV> r;
	r[0] = a[0][0]*b[0] + a[0][1]*b[1] + a[0][2];
	r[1] = a[1][0]*b[0] + a[1][1]*b[1] + a[1][2];
	return r;
}


//-----------------------------------------------------------------------------
//	TAffine
/**	TAffine�s��(��],�g��,���s�ړ���\��)�s��.
	�T�v�́C\ref pgAffine �Q�ƁD	*/

template <class T>
class TAffine:public PTM::TMatrixBase<4,4,
	PTM::TMatrixDescCol< TAffine<T>, PTM::TMatrixRow<4,4,T>, 4,4,4,T> >{
public:
	typedef PTM::TMatrixDescCol< TAffine<T>, PTM::TMatrixRow<4,4,T>, 4,4,4,T> desc;
	typedef PTM::TMatrixBase<4,4,desc> base_type;
	/**	�p������Ȃ���{�I�ȃ����o�̒�`.
		@see ::DEF_MATRIX_BASIC_MEMBER	*/
	DEF_MATRIX_BASIC_MEMBER(TAffine);
	union{
		struct{
			T xx, xy, xz, xw;
			T yx, yy, yz, yw;
			T zx, zy, zz, zw;
			T px, py, pz, pw;
		};
		T data[4][4];
	};
	///	�v�f�̃A�N�Z�X
	element_type& item_impl(size_t i, size_t j){ return data[j][i]; }
	const element_type& item_impl(size_t i, size_t j) const { return data[j][i]; }

	/**@name	���x�N�g���ւ̃A�N�Z�X	*/
	//@{
	/// 
	TVec3<element_type>& Ex() { return (TVec3<element_type>&)this->col(0); }
	/// 
	const TVec3<element_type>& Ex() const { return (TVec3<element_type>&)this->col(0); }
	/// 
	TVec3<element_type>& Ey() { return (TVec3<element_type>&)this->col(1); }
	/// 
	const TVec3<element_type>& Ey() const { return (TVec3<element_type>&)this->col(1); }
	/// 
	TVec3<element_type>& Ez() { return (TVec3<element_type>&)this->col(2); }
	/// 
	const TVec3<element_type>& Ez() const { return (TVec3<element_type>&)this->col(2); }
	/// 
	TVec3<element_type>& Trn() { return (TVec3<element_type>&)this->col(3); }
	/// 
	const TVec3<element_type>& Trn() const { return (TVec3<element_type>&)this->col(3); }
	///	���s�ړ�����(Trn()�̕ʖ�)
	TVec3<element_type>& Pos() {return Trn();}
	///	���s�ړ�����(Trn()�̕ʖ�,const ��)
	const TVec3<element_type>& Pos() const {return Trn();}
	//@}

	/**@name	�v�f�ւ̃A�N�Z�X	*/
	//@{
	/// 
	element_type& ExX() {return Ex().X();}
	const element_type& ExX() const {return Ex().X();}
	/// 
	element_type& ExY() {return Ex().Y();}
	const element_type& ExY() const {return Ex().Y();}
	/// 
	element_type& ExZ() {return Ex().Z();}
	const element_type& ExZ() const {return Ex().Z();}
	/// 
	element_type& EyX() {return Ey().X();}
	const element_type& EyX() const {return Ey().X();}
	/// 
	element_type& EyY() {return Ey().Y();}
	const element_type& EyY() const {return Ey().Y();}
	/// 
	element_type& EyZ() {return Ey().Z();}
	const element_type& EyZ() const {return Ey().Z();}
	/// 
	element_type& EzX() {return Ez().X();}
	const element_type& EzX() const {return Ez().X();}
	/// 
	element_type& EzY() {return Ez().Y();}
	const element_type& EzY() const {return Ez().Y();}
	/// 
	element_type& EzZ() {return Ez().Z();}
	const element_type& EzZ() const {return Ez().Z();}
	/// 
	element_type& TrnX() {return Trn().X();}
	const element_type& TrnX() const {return Trn().X();}
	/// 
	element_type& TrnY() {return Trn().Y();}
	const element_type& TrnY() const {return Trn().Y();}
	/// 
	element_type& TrnZ() {return Trn().Z();}
	const element_type& TrnZ() const {return Trn().Z();}
	///	TrnX()�̕ʖ�
	element_type& PosX() {return TrnX();}
	const element_type& PosX() const {return TrnX();}
	///	TrnY()�̕ʖ�
	element_type& PosY() {return TrnY();}
	const element_type& PosY() const {return TrnY();}
	///	TrnZ()�̕ʖ�
	element_type& PosZ() {return TrnZ();}
	const element_type& PosZ() const {return TrnZ();}
	///
	element_type& ExW() {return this->item(3,0);}
	const element_type& ExW() const {return this->item(3,0);}
	///
	element_type& EyW() {return this->item(3,1);}
	const element_type& EyW() const {return this->item(3,1);}
	///
	element_type& EzW() {return this->item(3,2);}
	const element_type& EzW() const {return this->item(3,2);}
	///
	element_type& TrnW() {return this->item(3,3);}
	const element_type& TrnW() const {return this->item(3,3);}
	///
	element_type& PosW() {return this->item(3,3);}
	const element_type& PosW() const {return this->item(3,3);}
	//@}

	///	��]�g��ϊ����ւ̎Q�Ƃ�Ԃ�.
	PTM::TSubMatrixCol<3,3, desc>& Rot(){ return this->sub_matrix(PTM::TSubMatrixDim<0,0,3,3>()); }
	///	��]�g��ϊ����ւ̎Q�Ƃ�Ԃ� (const��).
	const PTM::TSubMatrixCol<3,3, desc>& Rot() const { return this->sub_matrix(PTM::TSubMatrixDim<0,0,3,3>()); }

	/**@name	�������ƍ\�z	*/
	///	�P�ʍs��
	static TAffine<T> Unit(){
		TAffine<T> y;
		PTM::init_unitize(y);
		return y;
	}
	///	���s�ړ�
	static TAffine<T> Trn(element_type px, element_type py, element_type pz){
		TAffine<T> y;
		y.Trn().X() = px;
		y.Trn().Y() = py;
		y.Trn().Z() = pz;
		return y;
	}
	///	x/y/z���܂���]
	static TAffine<T> Rot(element_type th, char axis)
	{
		TAffine<T> y;
#ifdef __BORLANDC__
        TMatrix3<T> r = y.Rot();
		PTM::init_rot(r, th, axis);
#else
		PTM::init_rot(y.Rot(), th, axis);
#endif
		return y;
	}
	/**	�C�ӎ��܂���]
\verbatim
		+																	   +
		|u^2+(1-u^2)cos(th)      uv(1-cos(th))-wsin(th)  wu(1-cos(th))+vsin(th)|
	R =	|uv(1-cos(th))+wsin(th)  v^2+(1-v^2)cos(th)      vw(1-cos(th))-usin(th)|
		|wu(1-cos(th))-vsin(th)  vw(1-cos(th))+usin(th)  w^2+(1-w^2)cos(th)    |
		+																	   +
\endverbatim
*/
	template <class BUF>
	static TAffine<T> Rot(element_type th, const PTM::TVectorBase<3, BUF>& axis)
	{
		TAffine<T> y;
		Matrix3f r;
		PTM::init_rot(r, th, axis);
		y.Rot() = r;
		return y;
	}
	/// �g��
	static TAffine<T> Scale(element_type sx, element_type sy, element_type sz){
		TAffine<T> y;
		y.item(0, 0) = sx; y.item(1, 1) = sy; y.item(2, 2) = sz;
		return y;
	}
	/**	OpenGL�̎ˉe�s��Ƃ��ď�����
		@param screen	�J�������猩���r���[�|�[�g�̒��S�̈ʒu
		@param size		�r���[�|�[�g�̃T�C�Y
		@param front	��O�̃N���b�s���O���ʂƃJ�����̋���
		@param back		���̃N���b�s���O���ʂƃJ�����̋���	*/
	template <class BUFS, class BUFZ>
	static TAffine<T> ProjectionGL(
		const PTM::TVectorBase<3, BUFS>& screen,
		const PTM::TVectorBase<2, BUFZ>& size,
		element_type front=1.0f, element_type back=10000.0f)
	{
		TAffine<T> y;
		PTM::init_projection_gl(y, screen, size, front, back);
		return y;
	}
	/**	Direct3D�̎ˉe�s��Ƃ��ď�����
		@param screen	�J�������猩���r���[�|�[�g�̒��S�̈ʒu
		@param size		�r���[�|�[�g�̃T�C�Y
		@param front	��O�̃N���b�s���O���ʂƃJ�����̋���
		@param back		���̃N���b�s���O���ʂƃJ�����̋���	*/
	template <class BUFS, class BUFZ>
	static TAffine<T> ProjectionD3D(const PTM::TVectorBase<3, BUFS>& screen,
	const PTM::TVectorBase<2, BUFZ>& size,
		element_type front=1.0f, element_type back=10000.0f)
	{
		TAffine<T> y;
		PTM::init_projection_d3d(y, screen, size, front, back);
		return y;
	}

	/** OpenGL�̒����ˉe�ϊ�
		@param vpSize	�r���[�|�[�g�̃T�C�Y
	 */
	template <class BUFS>
	static TAffine<T> OrthoGL(const PTM::TVectorBase<2,BUFS>& vpSize){
		TAffine<T> y;
		PTM::init_ortho_gl(y, vpSize);
		return y;
	}

	/* obsolete. �f�t�H���g�����ɕς��܂��� tazz
	template <class BUF>
	void LookAt(const PTM::TVectorBase<3, BUF>& targetPos)
	{
		PTM::init_look_at(*this, targetPos);
	}
	*/
	/**	�����s��Ƃ��ď�����
		@param targetPos �����_�̈ʒu
		@param upDir �������x�N�g��
		Affine�s��̈ʒu�͂��̂܂܂ɁCEz()��targetPos�֌����CEy()��upDir�֌����悤�ɌX����ݒ肷��D
		�ʒu��Trn() (Pos())�ŗ\�ߐݒ肷��D
		// NOTE: ���܂ł̎����ł�upDir���ʒu�x�N�g���Ƃ��Ď�������Ă��܂������COpenGL�̎d�l�ɍ��킹�邽��
		// �����x�N�g���Ƃ��܂����D�����̃R�[�h�͉e�����󂯂�\��������܂��D tazz 09/05/29
	*/
	template <class BUF>
	void LookAt(const PTM::TVectorBase<3, BUF>& targetPos,
		const TVec3<typename BUF::element_type>& upDir = TVec3<typename BUF::element_type>(0,1,0))
	{
		PTM::init_look_at(*this, targetPos, upDir);
	}
	/* obsolete
	template <class BUF>
	void LookAtGL(const PTM::TVectorBase<3, BUF>& targetPos)
	{
		PTM::init_look_at_gl(*this, targetPos);
	}
	*/
	/** �����s��Ƃ��ď�����
		@param targetPos �����_�̈ʒu
		@param upDir �������x�N�g��
		Affine�s��̈ʒu�͂��̂܂܂ɁC-Ez()��targetPos�֌����CEy()��upDir�֌����悤�ɌX����ݒ肷��D
		�ʒu��Trn() (Pos())�ŗ\�ߐݒ肷��D
		// NOTE: ��Ɠ��l�̒��ӁD
	 */
	template <class BUF>
	void LookAtGL(const PTM::TVectorBase<3, BUF>& targetPos,
		const TVec3<typename BUF::element_type>& upDir = TVec3<typename BUF::element_type>(0,1,0))
	{
		PTM::init_look_at_gl(*this, targetPos, upDir);
	}
	/// ���K������
	void Orthonormalization(){
		this->Ex() /= this->Ex().norm();
		this->Ey() /= this->Ey().norm();
		this->Ez() /= this->Ez().norm();
		for (int i=0; i<10; ++i) {
			Vec3f u = PTM::cross(this->Ey(), this->Ez());
			Vec3f v = PTM::cross(this->Ez(), this->Ex());
			Vec3f w = PTM::cross(this->Ex(), this->Ey());
			this->Ex() = (this->Ex()+u)/2.0f;
			this->Ey() = (this->Ey()+v)/2.0f;
			this->Ez() = (this->Ez()+w)/2.0f;
			this->Ex() /= this->Ex().norm();
			this->Ey() /= this->Ey().norm();
			this->Ez() /= this->Ez().norm();
			float r = 0.0f;
			r += PTM::dot(this->Ex(), this->Ey())*PTM::dot(this->Ex(), this->Ey());
			r += PTM::dot(this->Ey(), this->Ez())*PTM::dot(this->Ey(), this->Ez());
			r += PTM::dot(this->Ez(), this->Ex())*PTM::dot(this->Ez(), this->Ex());
			if (r < 0.000001) {
				break;
			}
		}
	}
	
	///�R���X�g���N�^
	void set_default(){PTM::init_unitize(*this);}
};

/**	TAffine�̃R���X�g���N�^��`�̃}�N���D�R���X�g���N�^�͌p������Ȃ����C
	�����Define���邱�ƂŁC�R���X�g���N�^�������p�����Ƃ��ł���D	*/
#define DEF_TAFFINE_CONSTRUCTORS(TAffine)												\
	TAffine(){*this=Unit();}															\
	TAffine(element_type px, element_type py, element_type pz){*this=Trn(px, py, pz);}	\
	template <class BUFX, class BUFY>													\
	TAffine(const PTM::TVectorBase<3, BUFX>& exi,								\
			const PTM::TVectorBase<3, BUFY>& eyi){								\
			PTM::init_direct(Rot(), exi, eyi, 'x');										\
			item(3, 0) = 0; item(3, 1) = 0; item(3, 2) = 0; item(3, 3) = 1;				\
			item(0, 3) = 0; item(1, 3) = 0; item(2, 3) = 0;								\
	}																					\
	template <class BUFX, class BUFY, class BUFP>										\
	TAffine(	const PTM::TVectorBase<3, BUFX>& exi,							\
			const PTM::TVectorBase<3, BUFY>& eyi,								\
			const PTM::TVectorBase<3, BUFP>& posi){								\
			PTM::init_direct(Rot(), exi, eyi, 'x');										\
			item(3, 0) = 0; item(3, 1) = 0; item(3, 2) = 0; item(3, 3) = 1;				\
			item(0, 3) = posi.X(); item(1, 3) = posi.Y(); item(2, 3) = posi.Z();		\
	}                                                                                   \
	template <class BUFA, class BUFB>													\
	TAffine(	const PTM::TVectorBase<3, BUFA>& a,								\
			const PTM::TVectorBase<3, BUFB>& b,									\
			char axis){																	\
			PTM::init_direct(Rot(), exi, eyi, axis);									\
			item(3, 0) = 0; item(3, 1) = 0; item(3, 2) = 0; item(3, 3) = 1;				\
			item(0, 3) = 0; item(1, 3) = 0; item(2, 3) = 0;								\
	}																					\
	template <class BUFA, class BUFB, class BUFP>										\
	TAffine(	const PTM::TVectorBase<3, BUFA>& a,							    \
			const PTM::TVectorBase<3, BUFB>& b,									\
			char axis, const PTM::TVectorBase<3, BUFP>& posi){					\
			PTM::init_direct(Rot(), exi, eyi, axis);									\
			item(3, 0) = 0; item(3, 1) = 0; item(3, 2) = 0; item(3, 3) = 1;				\
			item(0, 3) = posi.X(); item(1, 3) = posi.Y(); item(2, 3) = posi.Z();		\
			}

/**	TAffine2�̃R���X�g���N�^��`�̃}�N���D�R���X�g���N�^�͌p������Ȃ����C
	�����Define���邱�ƂŁC�R���X�g���N�^�������p�����Ƃ��ł���D	*/
#define DEF_TAFFINE_CONSTRUCTORS2(TAffine)												\
	TAffine(element_type th, char axis) {												\
		*this=Rot(th, axis);															\
	}																					\
	template <class BUF>																\
	TAffine(element_type th, char axis,													\
		const PTM::TVectorBase<3, BUF>& posi) {									\
		*this = Rot(th, axis); Pos() = posi; }											\
	template <class BUFA>																\
	TAffine(element_type th, const PTM::TVectorBase<3, BUFA>& axis){			\
		*this = Rot(th, axis.unit());			                                        \
	}                                                                                   \
	template <class BUFA, class BUFP>													\
	TAffine(element_type th, const PTM::TVectorBase<3, BUFA>& axis , const PTM::TVectorBase<3, BUFP>& posi){	\
		*this = Rot(th, axis.unit()); Pos() = posi;										\
	}																					\
	template <class BUFS, class BUFZ>													\
	TAffine(const PTM::TVectorBase<3, BUFS>& screen, const PTM::TVectorBase<2, BUFZ>& size, element_type front=1.0f, element_type back=10000.0f){	\
		*this = ProjectionGL(screen, size, front, back);								\
	}                                                                                   \
	template <class BUF, class BUFV>													\
	TAffine(const PTM::TMatrixOp<3, 3, BUF>& m, const PTM::TVectorBase<3, BUFV> posi){\
		Rot() = m; Pos() = posi; ExW() = 0; EyW() = 0; EzW() = 0; PosW() = 1;			\
	}


#ifdef _WIN32
 #pragma warning (disable: 4700)
#endif

	
	
	///	TAffine�ƃx�N�g���̊|���Z
template <class TD, class TV>
TVec3<TV> operator * (const PTM::TMatrixBase<4,4, TD>& a, const TVec3<TV>& b){
	TVec3<TV> r;
	r[0] = a[0][0]*b[0] + a[0][1]*b[1] + a[0][2]*b[2] + a[0][3];
	r[1] = a[1][0]*b[0] + a[1][1]*b[1] + a[1][2]*b[2] + a[1][3];
	r[2] = a[2][0]*b[0] + a[2][1]*b[1] + a[2][2]*b[2] + a[2][3];
	return r;
}
#ifdef _WIN32
 #pragma warning (default: 4700)
#endif

///	float��2�����A�t�B���s��.
typedef TAffine2<float> Affine2f;
///	double��2�����A�t�B���s��.
typedef TAffine2<double> Affine2d;
///	float��3�����A�t�B���s��.
typedef TAffine<float> Affinef;
///	double��3�����A�t�B���s��.
typedef TAffine<double> Affined;
//@}
//@}


#ifndef PTM_PACK	//	�P�̂Ŏg�p����ꍇ�́Cnamespace �ɓ���Ȃ�
}	//	namespace Spr
#endif

#endif
