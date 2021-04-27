/*
 *  Copyright (c) 2003-2008, Shoichi Hasegawa and Springhead development team 
 *  All rights reserved.
 *  This software is free software. You can freely use, distribute and modify this 
 *  software. Please deal with this software under one of the following licenses: 
 *  This license itself, Boost Software License, The MIT License, The BSD License.   
 */
#ifndef BASE_DEBUG_H
#define BASE_DEBUG_H
#include "Env.h"
#include <stdarg.h>
#include <iostream>
#include <fstream>
#include <sstream>

/**	@file BaseDebug.h �f�o�b�O�p���[�e�B���e�B�[�N���X�E�֐��̒�`.	*/

#ifndef DBG_NEW
# ifdef _DEBUG
#  include <crtdbg.h>
#  define _CRTDBG_MAP_ALLOC
///	Debug�p new ���[�N���ɍs�ԍ���\��
#  define DBG_NEW  ::new(_NORMAL_BLOCK, __FILE__, __LINE__)
# else
#  define DBG_NEW new
# endif
#endif


/**	�f�o�b�O�p printf �֐�.
	@verbatim
	DPF("���b�Z�[�W:%s", msg);@endverbatim
	�̗l�Ɏg���D							*/
#define DPF	Spr::DebugPrintf::GetInstance()->FileLine(__FILE__, __LINE__)
/**	�f�o�b�O�p �o�̓X�g���[��.
	@verbatim
	DSTR << "���b�Z�[�W:" << msg;@endverbatim
	�̗l�Ɏg���D							*/
#define DSTR (Spr::DebugPrintf::GetInstance()->Stream())

#include <assert.h>

namespace Spr {

class SPR_DLL DebugPrintf{
public:
	static DebugPrintf* FASTCALL GetInstance();
	struct SPR_DLL PrintfFunc{
		const char* file;
		int line;
		PrintfFunc(const char* f, int l):file(f), line(l){}
		int SPR_CDECL operator() (const char*, ...);
	};
	PrintfFunc FileLine(const char* f=0, int l=-1){
		return PrintfFunc(f, l);
	}
	std::ostream& Stream();
	static void Set(void (*out)(const char*));
};

#if 0	//	�v���O�����̓�����ڍׂɕ񍐂�����Ȃ� 1
 #define TRACEALL DebugPrintf
#else
 #define TRACEALL (void*)
#endif

#if defined(_DEBUG) && !defined(NO_DEBUG_EVAL)
 #define DEBUG_EVAL(x) x
#else
 #define DEBUG_EVAL(x)
#endif


/**	�f�o�b�O�p CSV�o�̓X�g���[��.
	@verbatim
	CSVout << "���b�Z�[�W:" << msg;@endverbatim
	�̗l�Ɏg���D							*/
#define CSVOUT (Spr::DebugCSV::GetInstance()->Stream())
/**	�f�o�b�O�p CSV�o�̓X�g���[��.
	���݊J���Ă���t�@�C����close����D
	�ۑ�����ɂ͕K���ĂԕK�v������B
	*/
#define CSVCLOSE (Spr::DebugCSV::GetInstance()->Close())

class SPR_DLL DebugCSV{
public:
	static DebugCSV* instance;
	std::ofstream fout;

	static DebugCSV* FASTCALL GetInstance();
	std::ostream& Stream();
	void Set(void (*out)(const char*));
	static void defcsvOutFunc(const char* str);
	void Close();
	std::string FileNameSearch();							//�t�H���_����CSV�t�@�C�����T�[�`���V�����t�@�C�����𐶐�
};

}	//	namespace Spr

#endif
