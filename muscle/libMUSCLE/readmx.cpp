#include "libMUSCLE/muscle.h"
#include "libMUSCLE/textfile.h"

namespace muscle {

#define TRACE	0

const int MAX_LINE = 4096;
const int MAX_HEADINGS = 32;
static TLS<char[MAX_HEADINGS]> Heading;
static TLS<unsigned> HeadingCount(0);
static TLS<SCOREMATRIX> Mx;

static void LogMx()
	{
	Log("Matrix\n");
	Log("     ");
	for (int i = 0; i < 20; ++i)
		Log("    %c", LetterToChar(i));
	Log("\n");

	for (int i = 0; i < 20; ++i)
		{
		Log("%c    ", LetterToChar(i));
		for (int j = 0; j < 20; ++j)
			Log("%5.1f", Mx.get()[i][j]);
		Log("\n");
		}
	Log("\n");
	}

static unsigned MxCharToLetter(char c)
	{
	for (unsigned Letter = 0; Letter < HeadingCount.get(); ++Letter)
		if (Heading.get()[Letter] == c)
			return Letter;
	Quit("Letter '%c' has no heading", c);
	return 0;
	}

DYN_PTR_SCOREMATRIX ReadMx(TextFile &File)
	{
// Find column headers
	char Line[MAX_LINE];
	for (;;)
		{
		bool EndOfFile = File.GetLine(Line, sizeof(Line));
		if (EndOfFile)
			Quit("Premature EOF in matrix file");

		if (Line[0] == '#')
			continue;
		else if (Line[0] == ' ')
			break;
		else
			Quit("Invalid line in matrix file: '%s'", Line);
		}

// Read column headers
	HeadingCount.get() = 0;
	for (char *p = Line; *p; ++p)
		{
		char c = *p;
		if (!isspace(c))
			Heading.get()[HeadingCount.get()++] = c;
		}

	if (HeadingCount.get() > 0 && Heading.get()[HeadingCount.get()-1] == '*')
		--HeadingCount.get();

// AED 22/12/2006: don't force a 20 char matrix since nt matrices will have only 4
//	if (HeadingCount.get() < 20)
//		Quit("Error in matrix file: < 20 headers, line='%s'", Line);

#if TRACE
	{
	Log("ReadMx\n");
	Log("%d headings: ", HeadingCount.get());
	for (unsigned i = 0; i < HeadingCount.get(); ++i)
		Log("%c", Heading.get()[i]);
	Log("\n");
	}
#endif

// Zero out matrix
	for (int i = 0; i < MAX_ALPHA; ++i)
		for (int j = 0; j < MAX_ALPHA; ++j)
			Mx.get()[i][j] = 0.0;

// Read data lines
	for (unsigned RowIndex = 0; RowIndex < HeadingCount.get(); ++RowIndex)
		{
		bool EndOfFile = File.GetTrimLine(Line, sizeof(Line));
		if (EndOfFile)
			Quit("Premature EOF in matrix file");
#if	TRACE
		Log("Line=%s\n", Line);
#endif
		if (Line[0] == '#')
			continue;

		char c = Line[0];
#if	TRACE
		Log("Row char=%c\n", c);
#endif
		if (!IsResidueChar(c))
			continue;
		unsigned RowLetter = CharToLetter(c);
		if (RowLetter >= 20)
			continue;
#if	TRACE
		Log("Row letter = %u\n", RowLetter);
#endif

		char *p = Line + 1;
		char *maxp = p + strlen(Line);
		for (unsigned Col = 0; Col < HeadingCount.get() - 1; ++Col)
			{
			if (p >= maxp)
				Quit("Too few fields in line of matrix file: '%s'", Line);
			while (isspace(*p))
				++p;
			char *Value = p;
			while (!isspace(*p))
				++p;
			float v = (float) atof(Value);
			char HeaderChar = Heading.get()[Col];
			if (IsResidueChar(HeaderChar))
				{
				unsigned ColLetter = CharToLetter(HeaderChar);
				if (ColLetter >= 20)
					continue;
				Mx.get()[RowLetter][ColLetter] = v;
				}
			p += 1;
			}
		}

// Sanity check for symmetry
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < i; ++j)
			{
			if (Mx.get()[i][j] != Mx.get()[j][i])
				{
				Warning("Matrix is not symmetrical, %c->%c=%g, %c->%c=%g",
				  CharToLetter(i),
				  CharToLetter(j),
				  Mx.get()[i][j],
				  CharToLetter(j),
				  CharToLetter(i),
				  Mx.get()[j][i]);
				goto ExitLoop;
				}
			}
ExitLoop:;

	if (g_bVerbose.get())
		LogMx();

	SCOREMATRIX& sm = Mx.get();
	return &sm;
	}
} 
