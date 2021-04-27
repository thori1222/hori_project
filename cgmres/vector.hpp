
#pragma once

/// Springhead�̃x�N�g���E�s��N���X
#include <Base/TQuaternion.h>
//#include <initializer_list>

/// �x�N�g��
template<size_t n> struct vector_t_CSV_TEXT_GENERATOR;
template<size_t n>
struct vector_t : PTM::TVector<n, double>{
	template<class T> vector_t<n>& operator=(const T& v){ this->assign(v); return *this; }

	vector_t_CSV_TEXT_GENERATOR<n> ToCsv(){
		return vector_t_CSV_TEXT_GENERATOR<n>(*this);
	}

};

template <size_t n>
struct vector_t_CSV_TEXT_GENERATOR{
	friend vector_t<n>;

	vector_t<n>& trg;
	vector_t_CSV_TEXT_GENERATOR( vector_t<n>& t ):trg(t){};

private:	// noncopyable
	vector_t_CSV_TEXT_GENERATOR& operator=(const vector_t_CSV_TEXT_GENERATOR&);
	vector_t_CSV_TEXT_GENERATOR (const vector_t_CSV_TEXT_GENERATOR&);
};

template<size_t n>
std::ostream& operator<<(std::ostream &os, const vector_t_CSV_TEXT_GENERATOR<n>& trg)
{
	for( int i=0; i<n;i++){
		os << trg.trg[i];
		if(i != n-1)os << ",";
	}
	return os;
}


typedef Spr::Vec2d vec2_t;

