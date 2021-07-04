/*************************************************************************************************
CNFTools -- Copyright (c) 2015, Markus Iser, KIT - Karlsruhe Institute of Technology

--------------- Former Glucose Copyrights

 Glucose -- Copyright (c) 2009-2014, Gilles Audemard, Laurent Simon
                                CRIL - Univ. Artois, France
                                LRI  - Univ. Paris Sud, France (2009-2013)
                                Labri - Univ. Bordeaux, France

Glucose sources are based on MiniSat (see below MiniSat copyrights). Permissions and copyrights of
Glucose (sources until 2013, Glucose 3.0, single core) are exactly the same as Minisat on which it 
is based on. (see below).


--------------- Original Minisat Copyrights

Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*************************************************************************************************/

#ifndef SolverTypes_h
#define SolverTypes_h

#include <assert.h>
#include <stdint.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <iostream>

//=================================================================================================
// Variables, literals, lifted booleans:

struct Var {
	uint32_t id;

	inline Var(): id(0) { }
	inline Var(uint32_t id_): id(id_) { }

	inline operator int() const {
		return id;
	}

	inline Var& operator++ () { // prefix ++
		++id; return *this;
	}

	inline Var operator++ (int) { // postfix ++
		Var result(*this); ++id; return result;
	}
};

struct Lit {
	uint32_t x;

	inline Lit() : x(0) { }
	inline Lit(Var var) : Lit(var, false) { }
	inline Lit(Var var, bool sign) : x(2 * var + sign) { }

	inline operator int() const {
		return x;
	}

	inline bool sign() const {
		return x & 1;
	}

	inline Var var() const {
		return Var(x >> 1);
	}

	inline int toDimacs() {
		if (sign()) return (-1)*var(); else return var();
	}

	inline Lit operator ~ () const {
		Lit q;
		q.x = x ^ 1;
		return q;
	}

	inline Lit operator ^ (bool b) const {
		Lit q;
		q.x = x ^ (b ? 1 : 0);
		return q;
	}

	inline bool operator == (Lit p) const {
		return x == p.x;
	}

	inline bool operator != (Lit p) const {
		return x != p.x;
	}

	inline bool operator < (Lit p) const {
		return x < p.x;
	} // p and ~p will be adjacent
};

inline Var operator"" _V(unsigned long long n) {
    return Var((uint32_t)n);
}

inline Lit operator"" _L(unsigned long long n) {
    return Lit(Var((uint32_t)n));
}


const Var var_Undef(0);
const Lit lit_Undef(0, false);
const Lit lit_True(0, false);
const Lit lit_False(0, true);


//=================================================================================================
// Lifted booleans:

#define l_True  (lbool((uint8_t)0))
#define l_False (lbool((uint8_t)1))
#define l_Undef (lbool((uint8_t)2))

class lbool {
	uint8_t value;

public:
	explicit lbool(uint8_t v) : value(v) {}
	explicit lbool(bool x) : value(!x) {}	

	bool operator ==(lbool b) const {
		return (((b.value & 2) & (value & 2))
				| (!(b.value & 2) & (value == b.value))) != 0;
	}
	bool operator !=(lbool b) const {
		return !(*this == b);
	}
	lbool operator ^(bool b) const {
		return value == 2 ? l_Undef : lbool((uint8_t)(value ^ (uint8_t)b));
	}

	lbool operator &&(lbool b) const {
		uint8_t sel = (this->value << 1) | (b.value << 3);
		uint8_t v = (0xF7F755F4 >> sel) & 3;
		return lbool(v);
	}

	lbool operator ||(lbool b) const {
		uint8_t sel = (this->value << 1) | (b.value << 3);
		uint8_t v = (0xFCFCF400 >> sel) & 3;
		return lbool(v);
	}

	uint8_t operator |(lbool b) const {
		return this->value | b.value;
	}

	uint8_t operator &(lbool b) const {
		return this->value & b.value;
	}

	uint8_t operator |(uint8_t b) const {
		return this->value | b;
	}

	uint8_t operator &(uint8_t b) const {
		return this->value & b;
	}
};


typedef std::vector<Lit> Cl;
typedef std::vector<Cl*> For;

inline std::ostream& operator <<(std::ostream& stream, lbool const& value) {
	stream << (value == l_True ? '1' : (value == l_False ? '0' : 'X'));
    return stream;
}

inline std::ostream& operator <<(std::ostream& stream, Var const& var) {
	stream << var.id;
    return stream;
}

inline std::ostream& operator <<(std::ostream& stream, Lit const& lit) {
	if (lit.sign()) stream << "-";
	stream << lit.var();
    return stream;
}

inline std::ostream& operator <<(std::ostream& stream, Cl const& clause) {
    for (Lit lit : clause) {
        stream << lit << " ";
    }
    return stream;
}

inline std::ostream& operator <<(std::ostream& stream, For const& formula) {
    for (const Cl* clause : formula) {
        stream << *clause << std::endl;
    }
    return stream;
}

// add std::hash template specialization
namespace std {

template<>
struct hash<Lit> {
	std::size_t operator()(const Lit& key) const {
		Var hashedVar = key.var();
		if (key.sign() == 0) { 
			hashedVar = ~hashedVar;
		}
		return std::hash<uint32_t>()(hashedVar.id);
	}
};

template<>
struct hash<Var> {
	std::size_t operator()(const Var& key) const {
		return std::hash<uint32_t>()(key.id);
	}
};

}

#endif
