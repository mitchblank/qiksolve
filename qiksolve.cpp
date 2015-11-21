// qiksolve.cpp -- utility for finding solutions to qiktionary-style puzzles
//
//   Copyright (c) 2015 Mitchell Blank Jr
//   MIT licensed -- https://opensource.org/licenses/MIT
//
// On OS/X, compile with:
//   $ clang++ -std=c++11 -Wall -o qiksolve qiksolve.cpp
// I haven't tried other environments; some tweaking may be needed.
//
// Usage:
//   $ qiksolve %s [-t | -b] <word_file> [wordsize]
//
// "word_file" is a file of words, one per line.  The words that are not
// the correct length, have repeated characters, or have non-alphabetic
// characters are ignored.  On many UNIX machines, /usr/share/dict/words
// is a suitable input.
//
// If "wordsize" isn't given, defaults to 4.
//
// Flags:
//   -t -- Output a full solution tree for every word we know in
//         CSV format.  Fields are:
//           1. Vector of counts of matching letters from each step.  The
//              first move has an empty string
//           2. 'Y' if this is the last guess know how to make, 'N' otherwise
//           3. Suggestion to play.  If we have muliple anagrams, they're
//              given separated by vertical-bar characters
//           4. Bits of knowledge we will gain with this choice.  If this
//              is a final step this will be blank.  If this choice directly
//              leads to the answer it will be "-".
//           5. For non-final steps, a list of slash-separated counts of
//              how many possibilities we'll be left with if a given number
//              of matching characters match.
//
//   -b -- Batch mode builds the solution tree in memory, and then accepts
//         a series of words on standard input (one per line) and writes out
//         the steps we would play to solve it.
//
// If no flag is given, we'll go into interactive mode, where we present
// some suggestions to play.  Then you enter what you played and the number
// of characters that matched, and it will iteratively present a new list of
// suggestions.

#include <vector>
#include <set>
#include <map>
#include <string>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <cctype>
#include <cstring>

static int usage(char const * const argv0)
{
	fprintf(stderr, "%s: usage: %s [-t | -b] <word_file> [wordsize]\n", argv0, argv0);
	return 8;
}

#if defined(__GNUC__) || defined(__clang__)
#  define popcount32 	__builtin_popcount
#elif defined(_MSC_VER)
#  include <intrin.h>
#  define popcount32	__popcnt
#else
// Just a trivial implementation for other compilers
static unsigned popcount32(uint32_t const v)
{
	unsigned rv = 0;

	for (uint32_t b = 1; b != 0; b <<= 1)
		if ((v & b) != 0)
			rv++;
	return rv;
}
#endif

using namespace std;

// A "solution" can include more than one word in the case of anagrams
class Solution : public set<string> {
    public:
	Solution()
	{
	}
	void describe(FILE * const fout) const;
};

void Solution::describe(FILE * const fout) const
{
	assert(!empty());
	auto it = cbegin();
	fprintf(fout, "%s", it->c_str());
	while (++it != cend())
		fprintf(fout, "|%s", it->c_str());
}

// Since qiktionary does not allow repeated characters, we
// represent a word as a bitmask:
typedef uint32_t Word;

static Word ascii_to_bit(char const ch)
{
	Word b = 1;
	b <<= ((ch & 0x1F) - 1);
	return b;
}

static inline unsigned correct_letters(Word const a, Word const b)
{
	return popcount32(a & b);
}

typedef vector<Word> Words;

// When running in non-interactive mode, this tracks what results have been given
class AnswerSequence : public vector<unsigned> {
    public:
	AnswerSequence()
	{
	}
	void describe(FILE * const fout) const;
};

void AnswerSequence::describe(FILE * const fout) const
{
	if (!empty()) {
		auto it = cbegin();
		fprintf(fout, "%u", *it);
		while (++it != cend())
			fprintf(fout, ".%u", *it);
	}
}

class AbstractSolutionAcceptor;

// Class that holds all of the possible words
class Dictionary {
    public:
	Dictionary(FILE * const dictionary_file, unsigned const word_length);

	bool empty() const
	{
		return words_.empty();
	}

	// run in interactive mode:
	void interactive() const;

	// non-interactive mode:
	void printSolutions(FILE * const fout) const;

	// Read words on stdin, list how we would solve them on stdout
	void batch(FILE * const fin, FILE * const fout) const;

	unsigned wordLength() const
	{
		return word_length_;
	}

	// print out the word as a string, including any anagrmas
	void describeWord(FILE * const fout, Word const word) const
	{
		auto const it = words_.find(word);
		assert(it != words_.end());
		it->second.describe(fout);
	}

	// given a word, count how many anagrams we know
	size_t anagramCount(Word const word) const
	{
		auto const it = words_.find(word);
		assert(it != words_.end());
		return it->second.size();
	}
    private:
	// Non-interactive recursive solver
	void recursiveSolve(AbstractSolutionAcceptor * const acceptor, Words const& active_words, AnswerSequence * const runningAnswersp) const;
	// Interactive recursive suggestion engine
	void recursiveInteractive(Words const& active_words) const;
	// Print the best suggestions based on the entropy gained, given
	// what we know already
	void printInteractiveSuggestions(Words const& active_words, unsigned const show = 10) const;

	map<Word, Solution> words_;
	Words keys_;
	unsigned const word_length_;
};

// Build a bitmask "Word" one character at a time
class WordParser {
    public:
	WordParser()
		: accum_(0)
		, cnt_(0)
	{
	}
	// Take a single character
	bool parse(char const ch);
	// Get the result; 0 on failure
	Word result(unsigned const expected_length) const
	{
		return (cnt_ == expected_length) ? accum_ : 0;
	}
    private:
	Word accum_;
	unsigned cnt_;
};

bool WordParser::parse(char const ch)
{
	if (!isalpha(ch))
		return false;
	Word const bit = ascii_to_bit(ch);
	if ((accum_ & bit) != 0)
		return false;	// Words must have unique letters
	accum_ |= bit;
	cnt_++;
	return true;
}

// Given a string, turn it into a bitmask word
static Word string_to_bitmask(char const *p, unsigned const word_length)
{
	WordParser parser;

	while (*p != '\0')
		if (!parser.parse(*p++))
			return 0;
	return parser.result(word_length);
}

// Read a line from the file, turn it into uppercase
static bool read_uppercase_word(FILE * const fin, string *destp)
{
	destp->clear();
	for (;;) {
		int const ch = getc(fin);
		if (ch == EOF)
			return !destp->empty();
		if (ch == '\n')
			return true;
		destp->push_back(toupper(ch));
	}
}

Dictionary::Dictionary(FILE * const dictionary_file, unsigned const word_length)
	: word_length_(word_length)
{
	// Build words_[] map filled with all of the
	// valid words of "word_length" all-unique letters
	{
		string w;
		while (read_uppercase_word(dictionary_file, &w)) {
			Word const g = string_to_bitmask(w.c_str(), word_length);
			if (g != 0)
				words_[g].insert(w);
		}
	}
	if (!words_.empty()) {
		// Compute the keys_[] vector
		keys_.reserve(words_.size());
		size_t total_words = 0;
		for (auto it = words_.cbegin(); it != words_.cend(); ++it) {
			total_words += it->second.size();
			keys_.push_back(it->first);
		}
#if 0
		fprintf(stderr, "Seaching for %u (2^%lf) words (%u unique without anagrams)\n",
			static_cast<unsigned>(total_words),
			log(total_words) / M_LN2,
			static_cast<unsigned>(words_.size()));
#endif // 0
	}
}

// A vector of Word's along with the total number of choices
// they represent (including anagrams)
class WordsWithTotalChoices {
    public:
	WordsWithTotalChoices()
		: count_(0)
	{
	}
	void clear()
	{
		choices_.clear();
		count_ = 0;
	}
	size_t count() const
	{
		return count_;
	}
	// Add "guess" to the choices
	void update(Dictionary const& dictionary, Word const guess);

	Words const& choices() const
	{
		return choices_;
	}

    private:
	Words choices_;
	// count is different than choices.size() because we include anagrams
	size_t count_;
};

void WordsWithTotalChoices::update(Dictionary const& dictionary, Word const guess)
{
	choices_.push_back(guess);
	count_ += dictionary.anagramCount(guess);
}

// The wordlist split by the number of matching characters another word had
class DividedWordList : public vector<WordsWithTotalChoices> {
    public:
	explicit DividedWordList(unsigned const wordLength)
	{
		resize(wordLength + 1);
	}
	// Initialize the list, by splitting all of the words in "ws",
	// trying the word "guess"
	void build(Dictionary const& dictionary, Words const& ws, Word const guess);
	// Return how many bits of entropy are in this split (or a maximum
	// double if it perfectly splits between the various options)  We use
	// this to judge the "goodness" of a pick
	double entropy() const;
};

void DividedWordList::build(Dictionary const& dictionary, Words const& ws, Word const guess)
{
	assert(size() == dictionary.wordLength() + 1);
	assert(!ws.empty());
	for (auto it = begin(); it != end(); ++it)
		it->clear();
	{
		auto it = ws.cbegin();
		do {
			unsigned const c = correct_letters(*it, guess);
			assert(c <= dictionary.wordLength());
			(*this)[c].update(dictionary, *it);
		} while (++it != ws.cend());
	}
}

double DividedWordList::entropy() const
{
	size_t total_count = 0;
	bool perfect_split = true;
	for (auto it = cbegin(); it != cend(); ++it) {
		total_count += it->count();
		if (it->choices().size() > 1)
			perfect_split = false;
	}
	if (perfect_split) {
		// We can't do better than this -- by making this choice
		// we'll end up at the endgame
		return numeric_limits<double>::max();
	}
	double answer = 0.0;
	double const d_total = static_cast<double>(total_count);
	for (auto it = cbegin(); it != cend(); ++it)
		if (it->count() > 0) {
			double const prob = static_cast<double>(it->count()) / d_total;
			answer += prob * log(prob);
		}
	return -answer;
}

// Utility interface for Dictionary::recursiveSolve()
class AbstractSolutionAcceptor {
    public:
	virtual ~AbstractSolutionAcceptor();
	// Accept a non-terminal answer node, along with information
	// about the current entropy and information on how it is divided
	virtual void acceptBranch(AnswerSequence const& sequence, Word const word, double const entropy, DividedWordList const& divided) = 0;
	// Accept a terminal node -- we only have one more word that is
	// valid (although it can have anagrams)
	virtual void acceptSolution(AnswerSequence const& sequence, Word const word) = 0;
    protected:
	AbstractSolutionAcceptor()
	{
	}
};

AbstractSolutionAcceptor::~AbstractSolutionAcceptor()
{
	// nothing
}

void Dictionary::recursiveSolve(AbstractSolutionAcceptor * const acceptor, Words const& active_words, AnswerSequence * const runningAnswersp) const
{
	assert(!active_words.empty());
	if (active_words.size() == 1) {	// We found an end node
		acceptor->acceptSolution(*runningAnswersp, active_words.front());
	} else {
		DividedWordList dl(wordLength());
		double best_entropy = -1.0;
		Word best_word = 0;
		// Of all of the words in active_words, find the one that
		// will gives us the most even split (and thus the highest entropy)
		{
			auto it = keys_.cbegin();
			do {
				dl.build(*this, active_words, *it);
				double const entropy = dl.entropy();
				if (entropy > best_entropy) {
					best_entropy = entropy;
					best_word = *it;
				}
			} while (++it != keys_.cend());
		}
		assert(popcount32(best_word) == wordLength());
		// Now that we have found our best word, recompute the split
		dl.build(*this, active_words, best_word);
		assert(dl.entropy() == best_entropy);
		acceptor->acceptBranch(*runningAnswersp, best_word, best_entropy, dl);
		// Now recursively descend
		for (unsigned i = 0; i < wordLength() + 1; i++) {
			const WordsWithTotalChoices& wwtc = dl[i];
			if (wwtc.count() > 0) {
				runningAnswersp->push_back(i);
				recursiveSolve(acceptor, wwtc.choices(), runningAnswersp);
				runningAnswersp->pop_back();
			}
		}
	}
}

// Utility class for Dictionary::printSolutions().  For each step that
// recursiveSolve() finds, prints them to "fout"
class SolutionPrinter : public AbstractSolutionAcceptor {
    public:
	SolutionPrinter(Dictionary const& dict, FILE * const fout)
		: dict_(dict)
		, fout_(fout)
	{
	}
    private:
	void acceptBranch(AnswerSequence const& sequence, Word const word, double const entropy, DividedWordList const& divided) final;
	void acceptSolution(AnswerSequence const& sequence, Word const word) final;
	Dictionary const& dict_;
	FILE * const fout_;
};

void SolutionPrinter::acceptBranch(AnswerSequence const& sequence, Word const word, double const entropy, DividedWordList const& divided)
{
	sequence.describe(fout_);
	fputs(",N,", fout_);
	dict_.describeWord(fout_, word);
	// For debugging, print the entropy and the counts
	if (entropy == numeric_limits<double>::max())
		fputs(",-,", fout_);
	else
		fprintf(fout_, ",%lf,", entropy);
	{
		auto it = divided.cbegin();
		fprintf(fout_, "%u", static_cast<unsigned>(it->count()));
		while (++it != divided.cend())
			fprintf(fout_, "/%u", static_cast<unsigned>(it->count()));
	}
	putc('\n', fout_);
}

void SolutionPrinter::acceptSolution(AnswerSequence const& sequence, Word const word)
{
	sequence.describe(fout_);
	fputs(",Y,", fout_);
	dict_.describeWord(fout_, word);
	fputs(",,\n", fout_);
}

void Dictionary::printSolutions(FILE * const fout) const
{
	AnswerSequence seq;
	SolutionPrinter printer(*this, fout);
	recursiveSolve(&printer, keys_, &seq);
}

// Utility class for Dictionary::batch().  For each step that
// recursiveSolve() finds, fills them in a map for fast lookup.
class SolutionTree : public AbstractSolutionAcceptor {
    public:
	SolutionTree()
	{
	}
	Word operator[](string const& key) const
	{
		auto const it = map_.find(key);
		assert(it != map_.cend());
		return it->second;
	}
    private:
	void acceptBranch(AnswerSequence const& sequence, Word const word, double const entropy, DividedWordList const& divided) final;
	void acceptSolution(AnswerSequence const& sequence, Word const word) final;
	// map of "answer.answer."... to next suggestion
	map<string, Word> map_;
};

void SolutionTree::acceptBranch(AnswerSequence const& sequence, Word const word, double const, DividedWordList const&)
{
	acceptSolution(sequence, word);
}

static void append_unsigned_period(string * const destp, unsigned const u)
{
	char buf[100];
	snprintf(buf, sizeof(buf), "%u.", u);
	destp->append(buf);
}

void SolutionTree::acceptSolution(AnswerSequence const& sequence, Word const word)
{
	string key;
	for (auto it = sequence.cbegin(); it != sequence.cend(); ++it)
		append_unsigned_period(&key, *it);
	assert(map_.find(key) == map_.end());
	map_[key] = word;
}

void Dictionary::batch(FILE * const fin, FILE * const fout) const
{
	SolutionTree tree;
	{
		AnswerSequence seq;
		recursiveSolve(&tree, keys_, &seq);
	}
	string wstr;
	string key;
	while (read_uppercase_word(fin, &wstr)) {
		Word const target = string_to_bitmask(wstr.c_str(), wordLength());
		if (target == 0) {
			fprintf(stderr, "Bad word: %s\n", wstr.c_str());
		} else if (words_.find(target) == words_.end()) {
			fprintf(stderr, "Unknown  word: %s\n", wstr.c_str());
		} else {
			key.clear();
			printf("%s: ", wstr.c_str());
			for (;;) {
				Word const play = tree[key];
				assert(words_.find(play) != words_.end());
				if (play == target) {
					describeWord(fout, target);
					putc('\n', fout);
					break;
				}
				fputs(words_.find(play)->second.cbegin()->c_str(), fout);
				putc(' ', fout);
				append_unsigned_period(&key, correct_letters(target, play));
			}
		}
	}
}

// Parse a string in the form "<word> <matched>"
static bool parse_move(string const &cmd, unsigned const word_length, Word * const wordp, unsigned * const matchedp)
{
	WordParser parser;

	char const *p = cmd.c_str();
	for (;;) {
		if (*p == '\0')
			return false;
		if (*p == ' ' || *p == '\t')
			break;
		if (!parser.parse(*p++))
			return false;
	}
	*wordp = parser.result(word_length);
	if (*wordp == 0)
		return false;
	while (*p == ' ' || *p == '\t')
		p++;
	if (!isdigit(*p))
		return false;
	*matchedp = static_cast<unsigned>(atoi(p));
	return *matchedp <= word_length;
}

void Dictionary::recursiveInteractive(Words const& active_words) const
{
	if (active_words.empty()) {
		puts("No choices remain.\n");
	} else if (active_words.size() == 1) {
		fputs("Only choice remaining is: ", stdout);
		describeWord(stdout, active_words.front());
		puts("\n");
	} else {
		printf("%u choices remain", static_cast<unsigned>(active_words.size()));
		if (active_words.size() <= 5) {
			// If we're near the end, print a list
			fputs(": ", stdout);
			auto it = active_words.cbegin();
			describeWord(stdout, *it);
			do {
				fputs(", ", stdout);
				describeWord(stdout, *it);
			} while (++it != active_words.cend());
			putchar('\n');
		} else {
			puts(".");
		}
		putchar('\n');

		string cmd;
		for (;;) {
			printInteractiveSuggestions(active_words);
			fputs("Give your move in form \"<word> <result>\" like \"TAXI 2\"\n\n"
			      "> ", stdout);
			fflush(stdout);
			if (!read_uppercase_word(stdin, &cmd))
				break;	// EOF
			Word tried;
			unsigned matched;
			if (parse_move(cmd, wordLength(), &tried, &matched)) {
				DividedWordList dl(wordLength());
				dl.build(*this, active_words, tried);
				recursiveInteractive(dl[matched].choices());
				break;
			}
			puts("Invalid move.\n");
		}
	}
}

void Dictionary::printInteractiveSuggestions(Words const& active_words, unsigned const show) const
{
	assert(!active_words.empty());
	assert(show > 0);
	vector<pair<double, Word>> suggestions;
	suggestions.reserve(active_words.size());
	puts("Suggestions:");
	{
		DividedWordList dl(wordLength());
		auto it = keys_.cbegin();
		do {
			dl.build(*this, active_words, *it);
			suggestions.push_back(make_pair(-dl.entropy(), *it));
		} while (++it != keys_.cend());
	}
	sort(suggestions.begin(), suggestions.end());
	{
		auto it = suggestions.cbegin();
		unsigned i = 0;
		do {
			printf("\t%2u: ", ++i);
			describeWord(stdout, it->second);
			putchar('\n');
		} while (++it != suggestions.cend() && i < show);
	}
	putchar('\n');
}

void Dictionary::interactive() const
{
	recursiveInteractive(keys_);
}

int main(int argn, char const * const *argv)
{
	if (argn < 2)
		return usage(argv[0]);
	unsigned a = 1;
	enum {
		INTERACTIVE,
		SOLUTION_TREE,
		BATCH
	} mode = INTERACTIVE;
	if (0 == strcmp(argv[1], "-t")) {
		mode = SOLUTION_TREE;
		a = 2;
		argn--;
	} else if (0 == strcmp(argv[1], "-b")) {
		mode = BATCH;
		a = 2;
		argn--;
	}
	if (argn < 2 || argn > 3)
		return usage(argv[0]);
	FILE *fp = fopen(argv[a], "r");
	if (fp == NULL) {
		perror("open");
		return 4;
	}
	Dictionary dict(fp, (argn == 2) ? 4 : static_cast<unsigned>(atoi(argv[a + 1])));
	if (dict.empty()) {
		fprintf(stderr, "%s: No words found!\n", argv[0]);
		return 1;
	}
	switch (mode) {
	case INTERACTIVE:
		dict.interactive();
		break;
	case BATCH:
		dict.batch(stdin, stdout);
		break;
	default:
		assert(mode == SOLUTION_TREE);
		dict.printSolutions(stdout);
		break;
	}
	return 0;
}
