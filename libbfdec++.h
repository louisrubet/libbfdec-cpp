#pragma once
#include <complex>
#include <string>

extern "C" {
    #include "libbf.h"
}

typedef limb_t mp_prec_t;
// typedef bf_rnd_t mp_rnd_t;
typedef bf_flags_t mp_rnd_t;
typedef slimb_t mp_exp_t;

class Bfdec {
 private:
    bfdec_t mp;
    static bf_context_t bf_ctx;
    static void* bf_realloc(__attribute__((unused)) void*, void* ptr, size_t size) { return realloc(ptr, size); }

 public:
    // Get default rounding mode & precision
    inline static mp_rnd_t get_default_rnd() { return BF_RNDZ; }
    inline static mp_prec_t get_default_prec() { return BF_PREC_INF; }

    // Constructors && type conversions
    Bfdec() { bfdec_init(&bf_ctx, &mp); }
    Bfdec(const Bfdec& u) {
        bfdec_init(&bf_ctx, &mp);
        bfdec_set(&mp, &u.mp);
    }
    Bfdec(const double u, mp_prec_t prec = Bfdec::get_default_prec(), mp_rnd_t mode = Bfdec::get_default_rnd());
    Bfdec(const long double u, mp_prec_t prec = Bfdec::get_default_prec(), mp_rnd_t mode = Bfdec::get_default_rnd());
    Bfdec(const unsigned long long int u, mp_prec_t prec = Bfdec::get_default_prec(),
          mp_rnd_t mode = Bfdec::get_default_rnd());
    Bfdec(const long long int u, mp_prec_t prec = Bfdec::get_default_prec(), mp_rnd_t mode = Bfdec::get_default_rnd());
    Bfdec(const unsigned long int u, mp_prec_t prec = Bfdec::get_default_prec(),
          mp_rnd_t mode = Bfdec::get_default_rnd()) {
        bfdec_init(&bf_ctx, &mp);
        bfdec_set_ui(&mp, u);
    }
    Bfdec(const unsigned int u, mp_prec_t prec = Bfdec::get_default_prec(), mp_rnd_t mode = Bfdec::get_default_rnd());
    Bfdec(const long int u, mp_prec_t prec = Bfdec::get_default_prec(), mp_rnd_t mode = Bfdec::get_default_rnd()) {
        bfdec_init(&bf_ctx, &mp);
        bfdec_set_si(&mp, u);
    }
    Bfdec(const int u, mp_prec_t prec = Bfdec::get_default_prec(), mp_rnd_t mode = Bfdec::get_default_rnd()) {
        bfdec_init(&bf_ctx, &mp);
        bfdec_set_si(&mp, u);
    }

    // Construct Bfdec from bfdec_t structure.
    // shared = true allows to avoid deep copy, so that Bfdec and 'u' share the same data & pointers.
    Bfdec(const bfdec_t u, bool shared = false);

    Bfdec(const char* s, mp_prec_t prec = Bfdec::get_default_prec(), int base = 10,
          mp_rnd_t mode = Bfdec::get_default_rnd());
    Bfdec(const std::string& s, mp_prec_t prec = Bfdec::get_default_prec(), int base = 10,
          mp_rnd_t mode = Bfdec::get_default_rnd());

    ~Bfdec() { }

    static void clear_cache() { bf_clear_cache(&Bfdec::bf_ctx); }

#ifdef MPREAL_HAVE_MOVE_SUPPORT
    Bfdec& operator=(Bfdec&& v);
    Bfdec(Bfdec&& u);
#endif

    // Operations
    // =
    // +, -, *, /, ++, --, <<, >>
    // *=, +=, -=, /=,
    // <, >, ==, <=, >=

    // =
    Bfdec& operator=(const Bfdec& v) {
        bfdec_init(&bf_ctx, &mp);
        bfdec_set(&mp, &v.mp);
        return *this;
    }
    Bfdec& operator=(const long double v);
    Bfdec& operator=(const double v);
    Bfdec& operator=(const unsigned long int v);
    Bfdec& operator=(const unsigned long long int v);
    Bfdec& operator=(const long long int v);
    Bfdec& operator=(const unsigned int v);
    Bfdec& operator=(const long int v);
    Bfdec& operator=(const int v);
    Bfdec& operator=(const char* s);
    Bfdec& operator=(const std::string& s);
    template <typename real_t>
    Bfdec& operator=(const std::complex<real_t>& z);

    // +
    Bfdec& operator+=(const Bfdec& v) { return *this; }
    // Bfdec& operator+=(const mpf_t v);
    // Bfdec& operator+=(const mpz_t v);
    // Bfdec& operator+=(const mpq_t v);
    Bfdec& operator+=(const long double u);
    Bfdec& operator+=(const double u);
    Bfdec& operator+=(const unsigned long int u);
    Bfdec& operator+=(const unsigned int u);
    Bfdec& operator+=(const long int u);
    Bfdec& operator+=(const int u);

    Bfdec& operator+=(const long long int u);
    Bfdec& operator+=(const unsigned long long int u);
    Bfdec& operator-=(const long long int u);
    Bfdec& operator-=(const unsigned long long int u);
    Bfdec& operator*=(const long long int u);
    Bfdec& operator*=(const unsigned long long int u);
    Bfdec& operator/=(const long long int u);
    Bfdec& operator/=(const unsigned long long int u);

    const Bfdec operator+() const;
    Bfdec& operator++();
    const Bfdec operator++(int);

    // -
    Bfdec& operator-=(const Bfdec& v) {
        bfdec_sub(&mp, &mp, &v.mp, Bfdec::get_default_rnd(), Bfdec::get_default_prec());
        return *this;
    }
    // Bfdec& operator-=(const mpz_t v);
    // Bfdec& operator-=(const mpq_t v);
    Bfdec& operator-=(const long double u);
    Bfdec& operator-=(const double u);
    Bfdec& operator-=(const unsigned long int u);
    Bfdec& operator-=(const unsigned int u);
    Bfdec& operator-=(const long int u);
    Bfdec& operator-=(const int u);
    const Bfdec operator-() const {
        Bfdec r = *this;
        bfdec_neg(&r.mp);
        return r;
    }
    friend const Bfdec operator-(const unsigned long int b, const Bfdec& a);
    friend const Bfdec operator-(const unsigned int b, const Bfdec& a);
    friend const Bfdec operator-(const long int b, const Bfdec& a);
    friend const Bfdec operator-(const int b, const Bfdec& a);
    friend const Bfdec operator-(const double b, const Bfdec& a);
    Bfdec& operator--();
    const Bfdec operator--(int);

    // *
    Bfdec& operator*=(const Bfdec& v) {
        bfdec_mul(&mp, &mp, &v.mp, Bfdec::get_default_rnd(), Bfdec::get_default_prec());
        return *this;
    }
    // Bfdec& operator*=(const mpz_t v);
    // Bfdec& operator*=(const mpq_t v);
    Bfdec& operator*=(const long double v);
    Bfdec& operator*=(const double v);
    Bfdec& operator*=(const unsigned long int v);
    Bfdec& operator*=(const unsigned int v);
    Bfdec& operator*=(const long int v);
    Bfdec& operator*=(const int v);

    // /
    Bfdec& operator/=(const Bfdec& v) {
        bfdec_div(&mp, &mp, &v.mp, Bfdec::get_default_rnd(), Bfdec::get_default_prec());
        return *this;
    }
    // Bfdec& operator/=(const mpz_t v);
    // Bfdec& operator/=(const mpq_t v);
    Bfdec& operator/=(const long double v);
    Bfdec& operator/=(const double v);
    Bfdec& operator/=(const unsigned long int v);
    Bfdec& operator/=(const unsigned int v);
    Bfdec& operator/=(const long int v);
    Bfdec& operator/=(const int v);
    friend const Bfdec operator/(const unsigned long int b, const Bfdec& a);
    friend const Bfdec operator/(const unsigned int b, const Bfdec& a);
    friend const Bfdec operator/(const long int b, const Bfdec& a);
    friend const Bfdec operator/(const int b, const Bfdec& a);
    friend const Bfdec operator/(const double b, const Bfdec& a);

    //<<= Fast Multiplication by 2^u
    Bfdec& operator<<=(const unsigned long int u);
    Bfdec& operator<<=(const unsigned int u);
    Bfdec& operator<<=(const long int u);
    Bfdec& operator<<=(const int u);

    //>>= Fast Division by 2^u
    Bfdec& operator>>=(const unsigned long int u);
    Bfdec& operator>>=(const unsigned int u);
    Bfdec& operator>>=(const long int u);
    Bfdec& operator>>=(const int u);

    //
    inline bool gt0();
    inline bool ge0() { return mp.sign = 0; }
    inline bool lt0();
    inline bool le0();

    bool operator>(const Bfdec& v) const { return bfdec_cmp(&mp, &v.mp) > 0; }
    bool operator>=(const Bfdec& v) const { return bfdec_cmp(&mp, &v.mp) >= 0; }
    bool operator<(const Bfdec& v) const { return bfdec_cmp(&mp, &v.mp) < 0; }
    bool operator<=(const Bfdec& v) const { return bfdec_cmp(&mp, &v.mp) <= 0; }
    bool operator==(const Bfdec& v) const { return bfdec_cmp(&mp, &v.mp) == 0; }
    bool operator!=(const Bfdec& v) const { return bfdec_cmp(&mp, &v.mp) != 0; }

    // Type Conversion operators
    bfdec_t& toBfdec() { return mp; }
    const bfdec_t& toBfdec() const { return mp; }
    bool toBool() const;
    long toLong(mp_rnd_t mode = BF_RNDZ) const { return 0; }
    unsigned long toULong(mp_rnd_t mode = BF_RNDZ) const { return 0; }
    long long toLLong(mp_rnd_t mode = BF_RNDZ) const;
    unsigned long long toULLong(mp_rnd_t mode = BF_RNDZ) const;
    float toFloat(mp_rnd_t mode = BF_RNDN) const;
    double toDouble(mp_rnd_t mode = BF_RNDN) const;
    long double toLDouble(mp_rnd_t mode = BF_RNDN) const;

#if defined(MPREAL_HAVE_EXPLICIT_CONVERTERS)
    explicit operator bool() const { return toBool(); }
    explicit operator signed char() const { return (signed char)toLong(); }
    explicit operator unsigned char() const { return (unsigned char)toULong(); }
    explicit operator short() const { return (short)toLong(); }
    explicit operator unsigned short() const { return (unsigned short)toULong(); }
    explicit operator int() const { return (int)toLong(); }
    explicit operator unsigned int() const { return (unsigned int)toULong(); }
    explicit operator long() const { return toLong(); }
    explicit operator unsigned long() const { return toULong(); }
    explicit operator long long() const { return toLLong(); }
    explicit operator unsigned long long() const { return toULLong(); }
    explicit operator float() const { return toFloat(); }
    explicit operator double() const { return toDouble(); }
    explicit operator long double() const { return toLDouble(); }
#endif

    // Get raw pointers so that Bfdec can be directly used in raw mpfr_* functions
    // ::mpfr_ptr    mpfr_ptr();
    // ::mpfr_srcptr mpfr_ptr()    const;
    // ::mpfr_srcptr mpfr_srcptr() const;

    // Convert Bfdec to string with n significant digits in base b
    // n = -1 -> convert with the maximum available digits
    std::string toString(int n = -1, int b = 10, mp_rnd_t mode = Bfdec::get_default_rnd()) const;

    // #if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
    //     std::string toString(const std::string& format) const;
    // #endif

    std::ostream& output(std::ostream& os) const;

    // Math Functions
    friend const Bfdec sqr(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec sqrt(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec sqrt(const unsigned long int v, mp_rnd_t rnd_mode);
    friend const Bfdec cbrt(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec root(const Bfdec& v, unsigned long int k, mp_rnd_t rnd_mode);
    friend const Bfdec pow(const Bfdec& a, const Bfdec& b, mp_rnd_t rnd_mode);
    // friend const Bfdec pow (const Bfdec& a, const mpz_t b, mp_rnd_t rnd_mode);
    friend const Bfdec pow(const Bfdec& a, const unsigned long int b, mp_rnd_t rnd_mode);
    friend const Bfdec pow(const Bfdec& a, const long int b, mp_rnd_t rnd_mode);
    friend const Bfdec pow(const unsigned long int a, const Bfdec& b, mp_rnd_t rnd_mode);
    friend const Bfdec pow(const unsigned long int a, const unsigned long int b, mp_rnd_t rnd_mode);
    friend const Bfdec fabs(const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec abs(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec dim(const Bfdec& a, const Bfdec& b, mp_rnd_t rnd_mode);
    friend inline const Bfdec mul_2ui(const Bfdec& v, unsigned long int k, mp_rnd_t rnd_mode);
    friend inline const Bfdec mul_2si(const Bfdec& v, long int k, mp_rnd_t rnd_mode);
    friend inline const Bfdec div_2ui(const Bfdec& v, unsigned long int k, mp_rnd_t rnd_mode);
    friend inline const Bfdec div_2si(const Bfdec& v, long int k, mp_rnd_t rnd_mode);
    friend int cmpabs(const Bfdec& a, const Bfdec& b);

    friend const Bfdec log(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec log2(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec logb(const Bfdec& v, mp_rnd_t rnd_mode);
    friend mp_exp_t ilogb(const Bfdec& v);
    friend const Bfdec log10(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec exp(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec exp2(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec exp10(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec log1p(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec expm1(const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec nextpow2(const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec cos(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec sin(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec tan(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec sec(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec csc(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec cot(const Bfdec& v, mp_rnd_t rnd_mode);
    friend int sin_cos(Bfdec& s, Bfdec& c, const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec acos(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec asin(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec atan(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec atan2(const Bfdec& y, const Bfdec& x, mp_rnd_t rnd_mode);
    friend const Bfdec acot(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec asec(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec acsc(const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec cosh(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec sinh(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec tanh(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec sech(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec csch(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec coth(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec acosh(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec asinh(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec atanh(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec acoth(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec asech(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec acsch(const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec hypot(const Bfdec& x, const Bfdec& y, mp_rnd_t rnd_mode);

    friend const Bfdec fac_ui(unsigned long int v, mp_prec_t prec, mp_rnd_t rnd_mode);
    friend const Bfdec eint(const Bfdec& v, mp_rnd_t rnd_mode);

    friend const Bfdec gamma(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec tgamma(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec lngamma(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec lgamma(const Bfdec& v, int* signp, mp_rnd_t rnd_mode);
    friend const Bfdec zeta(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec erf(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec erfc(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec besselj0(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec besselj1(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec besseljn(long n, const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec bessely0(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec bessely1(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec besselyn(long n, const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec fma(const Bfdec& v1, const Bfdec& v2, const Bfdec& v3, mp_rnd_t rnd_mode);
    friend const Bfdec fms(const Bfdec& v1, const Bfdec& v2, const Bfdec& v3, mp_rnd_t rnd_mode);
    friend const Bfdec agm(const Bfdec& v1, const Bfdec& v2, mp_rnd_t rnd_mode);
    friend const Bfdec sum(const Bfdec tab[], const unsigned long int n, int& status, mp_rnd_t rnd_mode);
    friend int sgn(const Bfdec& v);

    // MPFR 2.4.0 Specifics
    // #if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))
    //     friend int          sinh_cosh   (Bfdec& s, Bfdec& c, const Bfdec& v, mp_rnd_t rnd_mode);
    //     friend const Bfdec li2         (const Bfdec& v,                       mp_rnd_t rnd_mode);
    //     friend const Bfdec fmod        (const Bfdec& x, const Bfdec& y,      mp_rnd_t rnd_mode);
    //     friend const Bfdec rec_sqrt    (const Bfdec& v,                       mp_rnd_t rnd_mode);
    //
    //     // MATLAB's semantic equivalents
    //     friend const Bfdec rem (const Bfdec& x, const Bfdec& y, mp_rnd_t rnd_mode); // Remainder after division
    //     friend const Bfdec mod (const Bfdec& x, const Bfdec& y, mp_rnd_t rnd_mode); // Modulus after division
    // #endif

    // #if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
    //     friend const Bfdec digamma (const Bfdec& v,        mp_rnd_t rnd_mode);
    //     friend const Bfdec ai      (const Bfdec& v,        mp_rnd_t rnd_mode);
    //     friend const Bfdec urandom (gmp_randstate_t& state, mp_rnd_t rnd_mode);     // use gmp_randinit_default() to
    //     init state, gmp_randclear() to clear
    // #endif

    // #if (MPFR_VERSION >= MPFR_VERSION_NUM(3,1,0))
    //     friend const Bfdec grandom (gmp_randstate_t& state, mp_rnd_t rnd_mode);     // use gmp_randinit_default() to
    //     init state, gmp_randclear() to clear friend const Bfdec grandom (unsigned int seed);
    // #endif

    // Uniformly distributed random number generation in [0,1] using
    // Mersenne-Twister algorithm by default.
    // Use parameter to setup seed, e.g.: random((unsigned)time(NULL))
    // Check urandom() for more precise control.
    friend const Bfdec random(unsigned int seed);

    // Splits Bfdec value into fractional and integer parts.
    // Returns fractional part and stores integer part in n.
    friend const Bfdec modf(const Bfdec& v, Bfdec& n);

    // Constants
    // don't forget to call mpfr_free_cache() for every thread where you are using const-functions
    friend const Bfdec const_log2(mp_prec_t prec, mp_rnd_t rnd_mode);
    friend const Bfdec const_pi(mp_prec_t prec, mp_rnd_t rnd_mode);
    friend const Bfdec const_euler(mp_prec_t prec, mp_rnd_t rnd_mode);
    friend const Bfdec const_catalan(mp_prec_t prec, mp_rnd_t rnd_mode);

    // returns +inf iff sign>=0 otherwise -inf
    friend const Bfdec const_infinity(int sign, mp_prec_t prec);

    // Output/ Input
    friend std::ostream& operator<<(std::ostream& os, const Bfdec& v);
    friend std::istream& operator>>(std::istream& is, Bfdec& v);

    // Integer Related Functions
    friend const Bfdec rint(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec ceil(const Bfdec& v);
    friend const Bfdec floor(const Bfdec& v);
    friend const Bfdec round(const Bfdec& v);
    friend long lround(const Bfdec& v);
    friend long long llround(const Bfdec& v);
    friend const Bfdec trunc(const Bfdec& v);
    friend const Bfdec rint_ceil(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec rint_floor(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec rint_round(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec rint_trunc(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec frac(const Bfdec& v, mp_rnd_t rnd_mode);
    friend const Bfdec remainder(const Bfdec& x, const Bfdec& y, mp_rnd_t rnd_mode);
    friend const Bfdec remquo(const Bfdec& x, const Bfdec& y, int* q, mp_rnd_t rnd_mode);

    // Miscellaneous Functions
    friend const Bfdec nexttoward(const Bfdec& x, const Bfdec& y);
    friend const Bfdec nextabove(const Bfdec& x);
    friend const Bfdec nextbelow(const Bfdec& x);

    // use gmp_randinit_default() to init state, gmp_randclear() to clear
    // friend const Bfdec urandomb (gmp_randstate_t& state);

    // MPFR < 2.4.2 Specifics
    // #if (MPFR_VERSION <= MPFR_VERSION_NUM(2,4,2))
    //     friend const Bfdec random2 (mp_size_t size, mp_exp_t exp);
    // #endif

    // Instance Checkers
    friend bool isnan(const Bfdec& v);
    friend bool isinf(const Bfdec& v);
    friend bool isfinite(const Bfdec& v);

    friend bool isnum(const Bfdec& v);
    friend bool iszero(const Bfdec& v);
    friend bool isint(const Bfdec& v);

    // #if (MPFR_VERSION >= MPFR_VERSION_NUM(3,0,0))
    //     friend bool isregular(const Bfdec& v);
    // #endif

    // Set/Get instance properties
    inline mp_prec_t get_prec() const;
    inline void set_prec(mp_prec_t prec, mp_rnd_t rnd_mode = get_default_rnd());  // Change precision with rounding mode

    // Aliases for get_prec(), set_prec() - needed for compatibility with std::complex<Bfdec> interface
    inline Bfdec& setPrecision(int Precision, mp_rnd_t RoundingMode = get_default_rnd());
    inline int getPrecision() const;

    // Set Bfdec to +/- inf, NaN, +/-0
    Bfdec& setInf(int Sign = +1);
    Bfdec& setNan();
    Bfdec& setZero(int Sign = +1);
    Bfdec& setSign(int Sign, mp_rnd_t RoundingMode = get_default_rnd());

    // Exponent
    mp_exp_t get_exp() const;
    int set_exp(mp_exp_t e);
    int check_range(int t, mp_rnd_t rnd_mode = get_default_rnd());
    int subnormalize(int t, mp_rnd_t rnd_mode = get_default_rnd());

    // Inexact conversion from float
    inline bool fits_in_bits(double x, int n);

    static limb_t bits2digits(mp_prec_t prec) {
        // const double LOG10_2 = 0.30102999566398119;
        // return int(std::floor( b * LOG10_2 ));
        return 16;
    }

    // Set/Get global properties
    static void set_default_prec(mp_prec_t prec) { } // TODO
    static void set_default_rnd(mp_rnd_t rnd_mode) { }

    static mp_exp_t get_emin(void);
    static mp_exp_t get_emax(void);
    static mp_exp_t get_emin_min(void);
    static mp_exp_t get_emin_max(void);
    static mp_exp_t get_emax_min(void);
    static mp_exp_t get_emax_max(void);
    static int set_emin(mp_exp_t exp);
    static int set_emax(mp_exp_t exp);

    // Efficient swapping of two Bfdec values - needed for std algorithms
    friend void swap(Bfdec& x, Bfdec& y);

    friend const Bfdec fmax(const Bfdec& x, const Bfdec& y, mp_rnd_t rnd_mode);
    friend const Bfdec fmin(const Bfdec& x, const Bfdec& y, mp_rnd_t rnd_mode);

 private:
    // Human friendly Debug Preview in Visual Studio.
    // Put one of these lines:
    //
    // mpfr::Bfdec=<DebugView>                              ; Show value only
    // mpfr::Bfdec=<DebugView>, <mp[0]._mpfr_prec,u>bits    ; Show value & precision
    //
    // at the beginning of
    // [Visual Studio Installation Folder]\Common7\Packages\Debugger\autoexp.dat
    // MPREAL_MSVC_DEBUGVIEW_DATA

    // "Smart" resources deallocation. Checks if instance initialized before deletion.
    // void clear(::mpfr_ptr);
};

typedef Bfdec mpreal;