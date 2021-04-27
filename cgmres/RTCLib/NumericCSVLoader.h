#pragma once

#ifndef RTCLIB_CSV_LOADER
#define RTCLIB_CSV_LOADER

//! ////////////////////////////////////////////////////////////////////////////////////
/// 
/// 自分のための，Realtime制御フレームワーク整備事業  2011 / Jun - 2011? / ??
///                                             Copyright @ 奥田裕之 2011 
/// 
/// NumericCSVLoader class : 簡便なCSV読み込み
///
/// 要件：C++, STL, Boost C++ library
/// 
/// ///////////////////////////////////////////////////////////////////////////////////////
///  
/// //---------------- 使い方 ----------------------
///
/// コンストラクタとしてデータファイル名を与えるだけ．
/// 二番目の引数には，ヘッダとして読み飛ばすべき行数を入れる．
/// 詳細なエラー処理が無いため，CSVとしての要件を守ること．
/// データ構造としては，行ベクトルを保持するベクトルになる．
///
/// //読み込み
///	RTCLib::NumericCSVLoader csv("demo_data001.csv",1);
///
/// //読み込んだデータへのアクセス
/// int i=0; // 行
/// int j=0; // 列
/// cout << csv[i][j];
/// 
/// -----------------[CSVファイルの要件]------------
/// ・ヘッダ行数が指定と正しいこと
/// ・コメント行などが無いこと
/// ・途中，データ行には文字が入っていないこと．
///
/// -----------------[その他注意]----------------------
/// Boost::Tokenizerを使っているので，大量のデータに対しては遅いかもしれない．
/// 
/// --------------------------------------------------------------------------------------

#include <iostream>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace RTCLib{

	// template<class SendObject, class RecvObject>
	class NumericCSVLoader
	{
	private:
		/// 読み込んだファイル名．フルパスに変換し格納するべし
		std::string m_filename;
		/// 読み込みが正常に行われているか
		bool m_valid;

		/// 列データ
		typedef std::vector<double> COL;
		typedef boost::shared_ptr<COL> pCOL;

		/// 保持データ
		std::vector< pCOL > m_data;

		/// ヘッダ文字列
		std::vector< std::string > m_headers;

		/// 縦横
		int m_rows;
		int m_cols;

	public:
		/// コンストラクタ
		/// 値の初期化のみ
		/// NumericCSVLoader();
		/// 読み込むべきファイル名つきコンストラクタ
		/// 初期化と同時に読み込む
		NumericCSVLoader();
		NumericCSVLoader(std::string fn, int header_lines);

		/// デストラクタ
		~NumericCSVLoader();

		/// 読み込みが正常に行われているか
		bool IsValid(){return m_valid;}

		/// ファイルを読み込む
		int Load(std::string fn, int header_lines);

		/// サイズ
		int Rows(){return m_rows;};
		int Cols(){return m_cols;};

		/// 表示
		void Print( std::ostream &os );

		/// 行要素へのアクセス．行超過に関しては，最終行が取得される．
		/// 行要素はvectorで定義されているので，NumericCSVLoader[i][j]でアクセス可能
		/// (iが行，jが列)
		const std::vector<double>& operator[](const int &n) const;

		/// ヘッダ文字列取得
		const std::string& Header(int n)const {return m_headers[n];}

		/// ヘッダ検索:列番号を返す．見つからなければ-1を返す．
		int GetColOf( const std::string &label );
		int GetColOf( const char *label );

	private:
		void Scan_Header( std::string headerstr );
		void Parse_Line( std::string &line, std::vector<double> &out );
		std::string trimstring(const std::string& s);
	};

}

#endif //RTCLIB_CSV_LOADER
