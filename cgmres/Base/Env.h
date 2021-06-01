/*
 *  Copyright (c) 2003-2008, Shoichi Hasegawa and Springhead development team 
 *  All rights reserved.
 *  This software is free software. You can freely use, distribute and modify this 
 *  software. Please deal with this software under one of the following licenses: 
 *  This license itself, Boost Software License, The MIT License, The BSD License.   
 */
#ifndef SPR_BASE_ENV_H
#define SPR_BASE_ENV_H

/**	@file Env.h
	�R���p�C���Ԃ̌݊�������邽�߂̃}�N����`.	*/
#define TYPENAME typename

#if defined __BORLANDC__
 #pragma warn -8026
 #pragma warn -8027
#endif


/**	DLL�̃G�N�X�|�[�g�̎w��
	DLL�ɃG�N�X�|�[�g���邽�߂ɂ́C
	class SPR_DLL DLLCLASS C{}; �̂悤�ɁC�N���X�錾�̑O�� DLLCLASS ������D
	DLL�����Ƃ��́CEXPORT_DLL ���}�N����`���Ă����D*/

#ifdef _MSC_VER
 #ifdef EXPORT_DLL
  #define SPR_DLL __declspec( dllexport )
  #pragma warning (disable: 4275 4251)
 #elif defined IMPORT_DLL
  #define SPR_DLL __declspec( dllimport )
 #else
  #define SPR_DLL
 #endif
#else
 #define SPR_DLL
#endif

/**	__cdecl
	Microsoft compiler �܂��� Borland C++ compiler �̏ꍇ�́A__cdecl���w��
 */
#if defined _MSC_VER || defined __BORLANDC__
 #define SPR_CDECL		__cdecl
 #define SPR_STDCALL 	__stdcall 
 #define FASTCALL		_fastcall
#else
 #define SPR_CDECL
 #define SPR_STDCALL 
 #define FASTCALL
#endif

/**	hdrstop
	�@Microsoft compiler�ABorland C++ compiler�A Intel C++ compiler �Ȃǂ̏ꍇ�́A 
	�@�v���R���p�C���w�b�_�[������s���B
 */
#ifndef __GNUC__
# define USE_HDRSTOP
#endif

/** ���l���Z�͈̓G���[�̃`�F�b�N
    VC���񋟂���֐��ł͈ꕔ�֐������قȂ�B
 */
#ifdef _MSC_VER
# include <float.h>
# define isnan  _isnan
# define finite _finite
#endif

//	for Visual C++ 's strange spec of stl.
#if defined _MSC_VER
 #if _MSC_VER <= 1300
  #undef min
  #undef max
  namespace std{
  template <class T>
  T min(T t1, T t2){
  	return t1 < t2 ? t1 : t2; 
  }
  template <class T>
  T max(T t1, T t2){
  	return t1 > t2 ? t1 : t2;
  }
  }
 #else
  #undef min
  #undef max
  #ifndef _MAX
   #define _MIN min
   #define _MAX max
  #endif
 #endif
#endif

//	x64�Ή�
#ifdef _WIN64
typedef unsigned __int64 ulong_ptr;
#else
typedef unsigned long ulong_ptr;
#endif

//	SWIG�p
#ifdef SWIGSPR
namespace PTM{
}
namespace Spr{
}
#endif

#endif
