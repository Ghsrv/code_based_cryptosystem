#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>
#include <cstring>


template < typename T >
class Matrix {
 private:
    size_t a, b;
    std::vector<T> matrix;
    const T& operator() (size_t x, size_t y) const {
        return matrix[x * b + y];
    }
    T& operator() (size_t x, size_t y) {
        return matrix[x * b + y];
    }
    Matrix(size_t a, size_t b): a(a), b(b) {
        (this->matrix).resize(a * b);
    }

 public:
    explicit Matrix(const std::vector <std::vector<T>>& matrix) {
        a = matrix.size();
        if (a != 0) {
            b = matrix[0].size();
        } else {
            b = 0;
        }
        (this->matrix).reserve(a * b);
        for (const auto& i : matrix) {
            for (const auto& j : i) {
                (this->matrix).push_back(j);
            }
        }
    }

    explicit Matrix(const std::vector<T>& matrix) {
        a = 1;
        b = matrix.size();
        (this->matrix).reserve(a * b);
        for (const auto& i : matrix) {
            (this->matrix).push_back(i);
        }
    }

    Matrix operator+ (const Matrix<T>& other) const {
        Matrix res = *this;
        for (size_t i = 0; i < matrix.size(); ++i) {
            res.matrix[i] += other.matrix[i];
        }
        return res;
    }

    Matrix& operator+= (const Matrix<T>& other) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            matrix[i] += other.matrix[i];
        }
        return *this;
    }

    template <typename t>
    Matrix& operator*= (const t& num) {
        for (size_t i = 0; i < matrix.size(); ++i) {
            matrix[i] *= num;
        }
        return *this;
    }

    template <typename t>
    Matrix operator* (const t& num) const {
        Matrix res = *this;
        for (size_t i = 0; i < matrix.size(); ++i) {
            res.matrix[i] *= num;
        }
        return res;
    }

    Matrix operator* (const Matrix& other) const {
        if (other.a != b) {
            throw std::runtime_error("Can't multiply matrices: wrong sizes " + std::to_string(b) + " != " + std::to_string(other.a));
        }
        Matrix res(a, other.b);
        for (size_t i = 0; i < a; ++i) {
            for (size_t j = 0; j < other.b; ++j) {
                for (size_t k = 0; k < b; ++k) {
                    res(i, j) += (*this)(i, k) * other(k, j);
                    // cout << i << ' ' << j << ' ' << res(i, j) << endl;
                }
            }
        }
        return res;
    }

    Matrix& operator*= (const Matrix& other) {
        *this = *this * other;
        return *this;
    }

    T& operator[] (size_t ind) {
        return matrix[ind];
    }

    Matrix transposed() const {
        Matrix t(b, a);
        for (size_t i = 0; i < a; ++i) {
            for (size_t j = 0; j < b; ++j) {
                t(j, i) = (*this)(i, j);
            }
        }
        return t;
    }

    Matrix& transpose() {
        *this = this->transposed();
        return *this;
    }

    auto begin() {
        return matrix.begin();
    }

    auto end() {
        return matrix.end();
    }

    auto begin() const {
        return matrix.cbegin();
    }

    auto end() const {
        return matrix.cend();
    }

    std::pair<size_t, size_t> size() const {
        return {a, b};
    }

    bool is_zero() const{
        return std::all_of(matrix.cbegin(), matrix.cend(), [](int i) { return i==0; });
    }

    friend std::ostream& operator<< (std::ostream& os, const Matrix<T>& a) {
        for (size_t i = 0; i < a.a; ++i) {
            for (size_t j = 0; j < a.b; ++j) {
                os << a(i, j);
                if (j != a.b - 1) os << '\t';
            }
            if (i != a.a - 1) os << '\n';
        }
        return os;
    }
};

class F2 {
 private:
    bool val;

 public:
    F2(): val(0) {}

    explicit F2(int val): val(val) {}

    operator bool() const {
        return this->val;
    }

    F2 operator*=(const F2& other) {
        this->val &= other.val;
        return *this;
    }

    F2 operator*(const F2& other) const {
        F2 res(this->val);
        return res *= other;
    }

    F2 operator+=(const F2& other) {
        this->val ^= other.val;
        return *this;
    }

    F2 operator+(const F2& other) const {
        F2 res(this->val);
        return res += other;
    }

    F2 operator=(const F2& other) {
        this->val = other.val;
        return *this;
    }

    F2 operator=(const int& val) {
        this->val = val;
        return *this;
    }

    friend std::ostream& operator<< (std::ostream& os, const F2& a) {
        os << a.val;
        return os;
    }

    friend std::istream& operator>> (std::istream& is, F2& a) {
        is >> a.val;
        return is;
    }
};


template<size_t s>
class Lfsr {
 private:
    std::bitset<s> bits;
    std::bitset<s> mask;

    int next() {
        int next = (bits & mask).count() % 2;
        bits <<= 1;
        bits.set(0, next);
        return next;
    }

 public:
    explicit Lfsr(std::bitset<s> bits, std::string mask): bits(bits), mask(mask) {}

    std::vector<F2> get() {
        std::vector<F2> res(s);
        for (auto &i : res) {
            i = this->next();
        }
        return res;
    }
};

F2 get_ins_bit(const std::vector<F2> v, size_t pos) {
    if (pos <= v.size() && v[pos - 1] == v[pos]) {
        return v[pos];
    }
    size_t cnt_1, cnt_0;
    for (int i = 1; i <= v.size() - pos; ++i) {
        if (i + pos == v.size() || v[pos + i] != v[pos]) {
            if (v[pos]) {
                cnt_1 = i;
            } else {
                cnt_0 = i;
            }
        }
    }
    for (int i = 1; i <= pos + 1; ++i) {
        if (i == pos + 1 || v[pos - i] != v[pos - 1]) {
            if (v[pos]) {
                cnt_1 = i;
            } else {
                cnt_0 = i;
            }
        }
    }
    if (cnt_1 > cnt_0) {
        return F2(1);
    } else if (cnt_1 < cnt_0) {
        return F2(0);
    }
    return F2(cnt_1 % 2);
}



int main(int argc, char *argv[]) {
    using namespace std;
    bool decode = 0;
    if (argc > 1 && strcmp(argv[1], "-d") == 0) {
        decode = 1;
    }
    ifstream pub;
    ifstream prk;
    pub.open("./keys/public_key");
    prk.open("./keys/private_key");
    bitset<90> rand1_state;
    bitset<128> rand2_state;
    bitset<128> rand3_state;
    for (int i = 0; i < 16; ++i) {
        char c;
        c = prk.get();
        for (int k = 0; k < 8; ++k) {
            if (c & 128) {
                rand2_state.set(i * 8 + k);
            }
            c <<= 1;
        }
    }
    for (int i = 0; i < 16; ++i) {
        char c;
        c = prk.get();
        for (int k = 0; k < 8; ++k) {
            if (c & 128) {
                rand3_state.set(i * 8 + k);
            }
            c <<= 1;
        }
    }
    for (int i = 0; i < 23; ++i) {
        char c;
        c = prk.get();
        for (int k = 0; (i * 8 + k < 90) && (k < 8) ; ++k) {
            if (c & 128) {
                rand1_state.set(i * 8 + k);
            }
            c <<= 1;
        }
    }

    prk.close();
    Lfsr<90> rand1(rand1_state, "100000000000000000000000000000000000000000000000000000000000000000000000000000000000010110");
    Lfsr<128> rand2(rand2_state, "10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000011");
    Lfsr<128> rand3(rand3_state, "10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001000011");

    vector<F2> r1 = rand1.get();
    vector<bool> insdel(15);
    vector<int> pos(15);
    int margin = 0;
    for (int i = 0; i < 15; ++i) {
        if (r1[i * 6]) {
            insdel[i] = 1;
            margin++;
        } else {
            margin--;
        }
        for (int j = 1; j < 6; ++j) {
            pos[i] <<= 1;
            if (r1[i * 6 + j]) {
                pos[i]++;
            }
        }
    }
    vector<vector<F2>> vG(256, vector<F2>(512));
    vector<vector<F2>> vH_1(512, vector<F2>(256));
    for (int i = 0; i < 256; ++i) {
        for (int j = 0; j < 64; ++j) {
            char c;
            c = pub.get();
            for (int k = 0; k < 8; ++k) {
                vG[i][j * 8 + k] = (c & 128) >> 7;
                c <<= 1;
            }
        }
    }
    for (int i = 0; i < 512; ++i) {
        for (int j = 0; j < 32; ++j) {
            char c;
            c = pub.get();
            for (int k = 0; k < 8; ++k) {
                vH_1[i][j * 8 + k] = ((c & 128) >> 7);
                c <<= 1;
            }
        }
    }
    pub.close();
    Matrix<F2> G(vG);
    Matrix<F2> G_T = G.transposed();
    Matrix<F2> H_1(vH_1);
    while (!cin.eof()) {
        vector<F2> r2 = rand2.get(), r3 = rand3.get();
        r2.insert(r2.end(), r3.begin(), r3.end());
        Matrix<F2> e = Matrix(r2) * H_1.transposed();
        char c = 0;
        int cnt = 0;

        if (!decode) {
            vector<vector<F2>> vm(1, vector<F2>(256));
            for (int i = 0; i < 32; ++i) {
                char c;
                c = cin.get();
                if (cin.eof()) {
                    break;
                }
                for (int k = 0; k < 8; ++k) {
                    vm[0][i * 8 + k] = (c & 128) >> 7;
                    c <<= 1;
                }
            }
            Matrix m(vm);
            Matrix res1 = m * G;
            res1 += e;
            vector<F2> ce(res1.begin(), res1.end());
            for (int i = 14; i >= 0; --i) {
                if (!insdel[i]) {
                    ce.erase(ce.begin() + i * (512 / 15) + pos[i]);
                } else {
                    ce.insert(ce.begin() + i * (512 / 15) + 1 + pos[i], get_ins_bit(ce, i * (512 / 15) + 1));
                }
            }
            for (auto i : ce) {
                c <<= 1;
                if (i) {
                    c++;
                }
                cnt++;
                if (cnt == 8) {
                    cout << c;
                    c = 0;
                    cnt = 0;
                }
            }
            if (cnt != 0) {
                c <<= 8 - cnt;
                cout << c;
            }
        } else {
            vector<F2> vc(512 + margin);
            for (int i = 0; i < (512 + margin + 7) / 8; ++i) {
                char c;
                c = cin.get();
                for (int k = 0; k < 8 && i * 8 + k < 512 + margin; ++k) {
                    vc[i * 8 + k] = (c & 128) >> 7;
                    c <<= 1;
                }
            }
            for (int i = 0; i < 15; ++i) {
                if (!insdel[i]) {
                    vc.insert(vc.begin() + i * (512 / 15) + pos[i], F2(0));
                } else {
                    vc.erase(vc.begin() + i * (512 / 15) + pos[i] + 1);
                }
            }
            vector<int> del_pos;
            for (int i = 0; i < 15; ++i) {
                if (!insdel[i]) {
                    del_pos.push_back((i * (512 / 15) + pos[i]));
                }
            }
            vector<vector<F2>> G_del_vec;
            for (auto i : del_pos) {
                G_del_vec.emplace_back(G_T.begin() + (256 * i), G_T.begin() + (256 * (i + 1)));
            }
            Matrix G_del(G_del_vec);
            Matrix ce(vc);
            ce += e;
            for (int i = 0; i < pos.size(); ++i) {
                if (!insdel[i]) {
                    G_del_vec.emplace_back(G_T.begin() + (256 * (i * (512 / 15) + pos[i])), G_T.begin() + (256 * (i * (512 / 15) + pos[i] + 1)));
                }
            }
            Matrix parity = ce * G_T;
            for (int i = 1; i < (1 << del_pos.size()); ++i) {
                vector<F2> v;
                for (int j = 0; j < del_pos.size(); ++j) {
                    v.push_back(F2((i >> j) & 1));
                }
                if ((Matrix(v) * G_del + parity).is_zero()) {
                    for (int j = 0; j < del_pos.size(); ++j) {
                        ce[del_pos[j]] = v[j] + ce[del_pos[j]];
                    }
                    break;
                }
            }
            ce *= H_1;
            for (auto i : ce) {
                c <<= 1;
                if (i) {
                    c++;
                }
                cnt++;
                if (cnt == 8) {
                    if (c == 0) {
                        break;
                    }
                    cout << c;
                    c = 0;
                    cnt = 0;
                }
            }
        }
        cin.peek();
    }
}
