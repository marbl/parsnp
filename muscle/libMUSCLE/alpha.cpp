#include "libMUSCLE/muscle.h"
#include <ctype.h>
#include "libMUSCLE/threadstorage.h"

namespace muscle {

/***
From Bioperl docs:
Extended DNA / RNA alphabet
------------------------------------------
Symbol       Meaning      Nucleic Acid
------------------------------------------
    A            A           Adenine
    C            C           Cytosine
    G            G           Guanine
    T            T           Thymine
    U            U           Uracil
    M          A or C
    R          A or G
    W          A or T
    S          C or G
    Y          C or T
    K          G or T
    V        A or C or G
    H        A or C or T
    D        A or G or T
    B        C or G or T
    X      G or A or T or C
    N      G or A or T or C

IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE:
         Cornish-Bowden (1985) Nucl. Acids Res. 13: 3021-3030.
***/
TLS<unsigned[MAX_CHAR]> g_CharToLetter;
TLS<unsigned[MAX_CHAR]> g_CharToLetterEx;

TLS<char[MAX_ALPHA]> g_LetterToChar;
TLS<char[MAX_ALPHA_EX]> g_LetterExToChar;

TLS<char[MAX_CHAR]> g_UnalignChar;
TLS<char[MAX_CHAR]> g_AlignChar;

TLS<bool[MAX_CHAR]> g_IsWildcardChar;
TLS<bool[MAX_CHAR]> g_IsResidueChar;

TLS<ALPHA> g_Alpha(ALPHA_Undefined);
TLS<unsigned> g_AlphaSize(0);

#define Res(c, Letter)												\
	{																\
	const unsigned char Upper = (unsigned char) toupper(c);			\
	const unsigned char Lower = (unsigned char) tolower(c);			\
	g_CharToLetter.get()[Upper] = Letter;									\
	g_CharToLetter.get()[Lower] = Letter;									\
	g_CharToLetterEx.get()[Upper] = Letter;								\
	g_CharToLetterEx.get()[Lower] = Letter;								\
	g_LetterToChar.get()[Letter] = Upper;									\
	g_LetterExToChar.get()[Letter] = Upper;								\
	g_IsResidueChar.get()[Upper] = true;									\
	g_IsResidueChar.get()[Lower] = true;									\
	g_AlignChar.get()[Upper] = Upper;										\
	g_AlignChar.get()[Lower] = Upper;										\
	g_UnalignChar.get()[Upper] = Lower;									\
	g_UnalignChar.get()[Lower] = Lower;									\
	}

#define Wild(c, Letter)												\
	{																\
	const unsigned char Upper = (unsigned char) toupper(c);			\
	const unsigned char Lower = (unsigned char) tolower(c);			\
	g_CharToLetterEx.get()[Upper] = Letter;								\
	g_CharToLetterEx.get()[Lower] = Letter;								\
	g_LetterExToChar.get()[Letter] = Upper;								\
	g_IsResidueChar.get()[Upper] = true;									\
	g_IsResidueChar.get()[Lower] = true;									\
	g_AlignChar.get()[Upper] = Upper;										\
	g_AlignChar.get()[Lower] = Upper;										\
	g_UnalignChar.get()[Upper] = Lower;									\
	g_UnalignChar.get()[Lower] = Lower;									\
	g_IsWildcardChar.get()[Lower] = true;									\
	g_IsWildcardChar.get()[Upper] = true;									\
	}

static unsigned GetAlphaSize(ALPHA Alpha)
	{
	switch (Alpha)
		{
	case ALPHA_Amino:
		return 20;

	case ALPHA_RNA:
	case ALPHA_DNA:
		return 4;
		}
	Quit("Invalid Alpha=%d", Alpha);
	return 0;
	}

static void InitArrays()
	{
	memset(g_CharToLetter.get(), 0xff, sizeof(g_CharToLetter.get()));
	memset(g_CharToLetterEx.get(), 0xff, sizeof(g_CharToLetterEx.get()));

	memset(g_LetterToChar.get(), '?', sizeof(g_LetterToChar.get()));
	memset(g_LetterExToChar.get(), '?', sizeof(g_LetterExToChar.get()));

	memset(g_AlignChar.get(), '?', sizeof(g_UnalignChar.get()));
	memset(g_UnalignChar.get(), '?', sizeof(g_UnalignChar.get()));

	memset(g_IsWildcardChar.get(), 0, sizeof(g_IsWildcardChar.get()));
	}

static void SetGapChar(char c)
	{
	unsigned char u = (unsigned char) c;

	g_CharToLetterEx.get()[u] = AX_GAP;
	g_LetterExToChar.get()[AX_GAP] = u;
	g_AlignChar.get()[u] = u;
	g_UnalignChar.get()[u] = u;
	}

static void SetAlphaDNA()
	{
	Res('A', NX_A)
	Res('C', NX_C)
	Res('G', NX_G)
	Res('T', NX_T)
	Wild('M', NX_M)
	Wild('R', NX_R)
	Wild('W', NX_W)
	Wild('S', NX_S)
	Wild('Y', NX_Y)
	Wild('K', NX_K)
	Wild('V', NX_V)
	Wild('H', NX_H)
	Wild('D', NX_D)
	Wild('B', NX_B)
	Wild('X', NX_X)
	Wild('N', NX_N)
	}

static void SetAlphaRNA()
	{
	Res('A', NX_A)
	Res('C', NX_C)
	Res('G', NX_G)
	Res('U', NX_U)
	Res('T', NX_T)
	Wild('M', NX_M)
	Wild('R', NX_R)
	Wild('W', NX_W)
	Wild('S', NX_S)
	Wild('Y', NX_Y)
	Wild('K', NX_K)
	Wild('V', NX_V)
	Wild('H', NX_H)
	Wild('D', NX_D)
	Wild('B', NX_B)
	Wild('X', NX_X)
	Wild('N', NX_N)
	}

static void SetAlphaAmino()
	{
	Res('A', AX_A)
	Res('C', AX_C)
	Res('D', AX_D)
	Res('E', AX_E)
	Res('F', AX_F)
	Res('G', AX_G)
	Res('H', AX_H)
	Res('I', AX_I)
	Res('K', AX_K)
	Res('L', AX_L)
	Res('M', AX_M)
	Res('N', AX_N)
	Res('P', AX_P)
	Res('Q', AX_Q)
	Res('R', AX_R)
	Res('S', AX_S)
	Res('T', AX_T)
	Res('V', AX_V)
	Res('W', AX_W)
	Res('Y', AX_Y)

	Wild('B', AX_B)
	Wild('X', AX_X)
	Wild('Z', AX_Z)
	}

void SetAlpha(ALPHA Alpha)
	{
	InitArrays();

	SetGapChar('.');
	SetGapChar('-');

	switch (Alpha)
		{
	case ALPHA_Amino:
		SetAlphaAmino();
		break;

	case ALPHA_DNA:
		SetAlphaDNA();

	case ALPHA_RNA:
		SetAlphaRNA();
		break;

	default:
		Quit("Invalid Alpha=%d", Alpha);
		}

	g_AlphaSize.get() = GetAlphaSize(Alpha);
	g_Alpha.get() = Alpha;

	if (g_bVerbose.get())
		Log("Alphabet %s\n", ALPHAToStr(g_Alpha.get()));
	}

char GetWildcardChar()
	{
	switch (g_Alpha.get())
		{
	case ALPHA_Amino:
		return 'X';

	case ALPHA_DNA:
	case ALPHA_RNA:
		return 'N';

	default:
		Quit("Invalid Alpha=%d", g_Alpha.get());
		}
	return '?';
	}

bool IsNucleo(char c)
	{
	return strchr("ACGTURYNacgturyn", c) != 0;
	}

bool IsDNA(char c)
	{
	return strchr("AGCTNagctn", c) != 0;
	}

bool IsRNA(char c)
	{
	return strchr("AGCUNagcun", c) != 0;
	}

static TLS<char[256]> InvalidLetters;
static TLS<int> InvalidLetterCount(0);

void ClearInvalidLetterWarning()
	{
	memset(InvalidLetters.get(), 0, 256);
	}

void InvalidLetterWarning(char c, char w)
	{
	InvalidLetters.get()[(unsigned char) c] = 1;
	++InvalidLetterCount.get();
	}

void ReportInvalidLetters()
	{
	if (0 == InvalidLetterCount.get())
		return;

	char Str[257];
	memset(Str, 0, 257);

	int n = 0;
	for (int i = 0; i < 256; ++i)
		{
		if (InvalidLetters.get()[i])
			Str[n++] = (char) i;
		}
	Warning("Assuming %s (see -seqtype option), invalid letters found: %s",
	  ALPHAToStr(g_Alpha.get()), Str);
	}
} 
