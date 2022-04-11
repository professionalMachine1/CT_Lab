#pragma once

#ifndef _NEXTGENMATRIX_H
#define  _NEXTGENMATRIX_H

#include <initializer_list>
#include <iostream>
#include <complex>
#include <math.h>

template <class ValType>
class Vector
{
protected:
	ValType* pVector;
	int size_;
public:
	//Конструкторы
	Vector(int _size = 1);
	Vector(const Vector& Vec);
	Vector(const ValType& Val, int _size);
	Vector(std::initializer_list<ValType> in_list);

	//Оператор присваивания
	Vector& operator= (const Vector& Vec);
	Vector& operator= (std::initializer_list<ValType> in_list);

	//Операции с векторами
	Vector  operator+(const Vector& Vec);
	Vector  operator-(const Vector& Vec);
	ValType  operator*(const Vector& Vec);
	Vector operator+=(const Vector& Vec);
	Vector operator-=(const Vector& Vec);

	//Операции с числами
	Vector  operator*(const ValType& Val);
	Vector  operator/(const ValType& Val);

	//Размер массива
	int size() const { return size_; }

	//Взятие значения по индексу
	//return (ind < 0) || (ind >= size) ? ValType() : pVector[i];
	ValType& operator[](int ind) const { return pVector[ind]; }

	//Логические операции
	bool operator==(const Vector& Vec) const;
	bool operator!=(const Vector& Vec) const;
	bool operator>(const Vector& Vec) const;
	bool operator<(const Vector& Vec) const;
	bool operator>=(const Vector& Vec) const;
	bool operator<=(const Vector& Vec) const;

	//Инкремент и декремент
	Vector& operator++();
	Vector operator++(int);
	Vector& operator--();
	Vector operator--(int);

	//Функции ввода-выводы
	friend std::istream& operator >>(std::istream& is, Vector& Vec)
	{
		for (int i = 0; i < Vec.size_; i++)
			is >> Vec.pVector[i];
		return is;
	}
	friend std::ostream& operator<< (std::ostream& os, const Vector& Vec)
	{
		for (int i = 0; i < Vec.size_; i++)
			os << Vec.pVector[i] << "   ";
		os << "\n";
		return os;
	}

	~Vector();

	//Также вне класса определено несколько полезных методов
	//Vector<ValType> operator*(const ValType& Val, const Vector<ValType>& Vec)
	//int MaxValue(const Vector<ValType>& Vec)
	//int MinValue(const Vector<ValType>& Vec)
	//Vector<ValType> abs(const Vector<ValType>& Vec)
	//ValType operator*(const ValType& Val1, const ValType& Val2)

};

/*

	Конструкторы

*/

template <class ValType>
Vector<ValType>::Vector(int _size)
{
	size_ = _size;
	pVector = new ValType[size_];
	for (int i = 0; i < size_; i++)
		pVector[i] = ValType();
}

template <class ValType>
Vector<ValType>::Vector(const ValType& Val, int _size)
{
	size_ = _size;
	pVector = new ValType[size_];
	for (int i = 0; i < size_; i++)
		pVector[i] = Val;
}

template <class ValType>
Vector<ValType>::Vector(const Vector& Vec)
{
	size_ = Vec.size_;
	pVector = new ValType[size_];
	for (int i = 0; i < size_; i++)
		pVector[i] = Vec.pVector[i];
}

template <typename ValType>
Vector<ValType>::Vector(std::initializer_list<ValType> in_list)
{
	size_ = in_list.size();
	pVector = new ValType[size_];

	auto it = in_list.begin();
	for (int i = 0; i < size_; i++)
		pVector[i] = *it, it++;
}

/*

	Оператор присваивания

*/

template <class ValType>
Vector<ValType>& Vector<ValType>::operator= (const Vector<ValType>& Vec)
{
	if (this != &Vec)
	{
		if (size_ != Vec.size_)
		{
			delete[]pVector;
			size_ = Vec.size_;
			pVector = new ValType[size_];
		}
		for (int i = 0; i < size_; i++)
			pVector[i] = Vec.pVector[i];
	}
	return *this;
}

template <class ValType>
Vector<ValType>& Vector<ValType>::operator= (std::initializer_list<ValType> in_list)
{
	if (size_ != in_list.size())
	{
		delete[]pVector;
		size_ = in_list.size();
		pVector = new ValType[size_];
	}

	auto it = in_list.begin();
	for (int i = 0; i < size_; i++)
		pVector[i] = *it, it++;
	return *this;
}

/*

	Операции с векторами

*/

template <class ValType>
Vector<ValType>  Vector<ValType>::operator+(const Vector<ValType>& Vec)
{
	//int n = size < Vec.size ? size : Vec.size;
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] + Vec.pVector[i];
	return res;
}

template <class ValType>
Vector<ValType>  Vector<ValType>::operator-(const Vector<ValType>& Vec)
{
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] - Vec.pVector[i];
	return res;
}

template <class ValType>
ValType  Vector<ValType>::operator*(const Vector<ValType>& Vec)
{
	ValType res = ValType();
	for (int i = 0; i < size_; i++)
		res += pVector[i] * Vec.pVector[i];
	return res;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator+=(const Vector<ValType>& Vec)
{
	Vector<ValType> res = *this;
	for (int i = 0; i < size_; i++)
		pVector[i] += Vec.pVector[i];
	return res;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator-=(const Vector<ValType>& Vec)
{
	Vector<ValType> res = *this;
	for (int i = 0; i < size_; i++)
		pVector[i] -= Vec.pVector[i];
	return res;
}

/*

	Операции с числами

*/

template <class ValType>
Vector<ValType>  Vector<ValType>::operator*(const ValType& Val)
{
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] * Val;
	return res;
}

template <class ValType>
Vector<ValType>  Vector<ValType>::operator/(const ValType& Val)
{
	Vector<ValType> res(size_);
	for (int i = 0; i < size_; i++)
		res.pVector[i] = pVector[i] / Val;
	return res;
}

//template <class ValType>
//ValType operator*(const ValType& Val1, const ValType& Val2)
//{
//	return Val1 * Val2;
//}

template <class ValType>
Vector<ValType> operator*(const ValType& Val, const Vector<ValType>& Vec)
{
	Vector<ValType> res(Vec.size());
	for (int i = 0; i < Vec.size(); i++)
		res[i] = Vec[i] * Val;
	return res;
}

/*

	Логические операции

*/

template <class ValType>
bool Vector<ValType>::operator==(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] != Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator!=(const Vector<ValType>& Vec) const
{
	return !(*this == Vec);
}

template <class ValType>
bool Vector<ValType>::operator>(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] <= Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator<(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] >= Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator>=(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] < Vec.pVector[i])
			return false;
	return true;
}

template <class ValType>
bool Vector<ValType>::operator<=(const Vector<ValType>& Vec) const
{
	if (size_ != Vec.size_)
		return false;
	for (int i = 0; i < size_; i++)
		if (pVector[i] > Vec.pVector[i])
			return false;
	return true;
}

/*

	Инкремент и декремент

*/

template <class ValType>
Vector<ValType>& Vector<ValType>::operator++()
{
	for (int i = 0; i < size_; i++)
		++pVector[i];
	return *this;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator++(int)
{
	Vector<ValType> res = *this;
	++* this;
	return res;
}

template <class ValType>
Vector<ValType>& Vector<ValType>::operator--()
{
	for (int i = 0; i < size_; i++)
		--pVector[i];
	return *this;
}

template <class ValType>
Vector<ValType> Vector<ValType>::operator--(int)
{
	Vector<ValType> res = *this;
	--* this;
	return res;
}

/*

	Полезные операции

*/

template <class ValType>
int MaxValue(const Vector<ValType>& Vec)
{
	ValType Max = Vec[0];
	int IndMax = 0;
	for (int i = 1; i < Vec.size(); i++)
		if (Vec[i] > Max)
		{
			Max = Vec[i];
			IndMax = i;
		}
	return IndMax;
}

template <class ValType>
int MinValue(const Vector<ValType>& Vec)
{
	ValType Min = Vec[0];
	int IndMin = 0;
	for (int i = 1; i < Vec.size(); i++)
		if (Vec[i] < Min)
		{
			Min = Vec[i];
			IndMin = i;
		}
	return IndMin;
}

template <class ValType>
Vector<ValType> abs(const Vector<ValType>& Vec)
{
	Vector<ValType> res(Vec.size());
	for (int i = 0; i < Vec.size(); i++)
		res[i] = fabs(Vec[i]);
	return res;
}

/*

	Деструктор

*/

template <class ValType>
Vector<ValType>::~Vector()
{
	delete[]pVector;
}

/*



	Класс матриц



*/

template <class ValType>
class Matrix : public Vector<Vector<ValType> >
{
private:
	int col_ind;
public:
	//Конструкторы 
	Matrix(int _size = 1);
	Matrix(const Matrix& mt);
	Matrix(int rows, int cells);
	Matrix(const Vector<Vector<ValType> >& mt);
	Matrix(const ValType& Val, int rows, int cells);
	Matrix(std::initializer_list<std::initializer_list<ValType>> in_list);

	//Оператор присваивания
	Matrix& operator= (const Matrix& mt);
	Matrix& operator= (std::initializer_list<std::initializer_list<ValType>> in_list);
	Matrix& operator= (const Vector<ValType>& Vec);
	Matrix& operator= (std::initializer_list<ValType> init_list);

	//Взятие по индексу
	ValType& operator()(int i, int j) { return this->pVector[i][j]; }
	Matrix& operator()(int j) { col_ind = j; return *this; }

	//Полезные операции
	Matrix t();
	Matrix i();
	Matrix& diag(const Vector<ValType>& diag_el);
	Matrix& diag(std::initializer_list<ValType> init_list);

	//Операции с матрицами
	Matrix operator*(const Matrix& mt);
	Matrix operator+(const Matrix& mt);
	Matrix operator-(const Matrix& mt);
	Matrix operator+=(const Matrix& mt);
	Matrix operator-=(const Matrix& mt);
	Matrix operator*=(const Matrix& mt);

	//Операции с числами
	Matrix  operator*(const ValType& Val);
	Matrix  operator/(const ValType& Val);

	//Логические операции
	bool operator==(const Matrix& mt) const;
	bool operator!=(const Matrix& mt) const;
	bool operator>(const Matrix& mt) const;
	bool operator<(const Matrix& mt) const;
	bool operator>=(const Matrix& mt) const;
	bool operator<=(const Matrix& mt) const;

	//Инкремент и декремент
	Matrix& operator++();
	Matrix operator++(int);
	Matrix& operator--();
	Matrix operator--(int);

	friend std::istream& operator >> (std::istream& is, Matrix& mt)
	{
		for (int i = 0; i < mt.size_; i++)
			is >> mt.pVector[i];
		return is;
	}
	friend std::ostream& operator << (std::ostream& os, const Matrix& mt)
	{
		os << std::endl;
		for (int i = 0; i < mt.size_; i++)
			os << mt.pVector[i];
		return os;
	}

	//Несколько полезных методов вне класса
	//Vector<ValType> operator*(const Vector<ValType>& Vec, const Matrix<ValType>& mt)
	//Vector<ValType> operator*(const Matrix<ValType>& mt, const Vector<ValType>& Vec)
	//Matrix<ValType> Single(int size_)
	//Matrix<ValType> pow(const Matrix<ValType>& mt, const int& degree)
	//Matrix<ValType> InverseMatrix(const Matrix<ValType>& mt)
	//Matrix<ValType> TransposedMatrix(const Matrix<ValType>& mt)
	//ValType Determinant(const Matrix<ValType>& mt)
	//Matrix<ValType> Minor(const Matrix<ValType>& mt, int ntrow, int ntcell)
	//Matrix<ValType> CellRowMultiply(const Vector<ValType>& cell, const Vector<ValType>& row)
	//Matrix<ValType> operator*(const ValType& Val, const Matrix<ValType>& mt)

};

/*

	Конструкторы

*/

template <class ValType>
Matrix<ValType>::Matrix(int _size) : Vector<Vector<ValType>>(_size)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = Vector<ValType>(Matrix::size_);
}

template <class ValType>
Matrix<ValType>::Matrix(const Matrix<ValType>& mt) : Vector<Vector<ValType>>(mt) { }

template <class ValType>
Matrix<ValType>::Matrix(int rows, int cells) : Vector<Vector<ValType>>(rows)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = Vector<ValType>(cells);
}

template <class ValType>
Matrix<ValType>::Matrix(const Vector<Vector<ValType> >& mt) : Vector<Vector<ValType>>(mt) { }


template <class ValType>
Matrix<ValType>::Matrix(const ValType& Val, int rows, int cells) : Vector<Vector<ValType>>(rows)
{
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = Vector<ValType>(Val, cells);
}

template <class ValType>
Matrix<ValType>::Matrix(std::initializer_list<std::initializer_list<ValType>> in_list)
{
	Matrix::size_ = in_list.size();
	Matrix::pVector = new Vector<ValType>[Matrix::size_];

	auto it = in_list.begin();
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = *it, it++;
}

/*

	Оператор присваивания

*/

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator= (const Matrix& mt)
{
	if (this != &mt)
	{
		if (Matrix::size_ != mt.size_)
		{
			delete[]Matrix::pVector;
			Matrix::size_ = mt.size_;
			Matrix::pVector = new Vector<ValType>[Matrix::size_];
		}
		for (int i = 0; i < Matrix::size_; i++)
			Matrix::pVector[i] = mt.pVector[i];
	}
	return *this;
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator= (std::initializer_list<std::initializer_list<ValType>> in_list)
{
	if (Matrix::size_ != in_list.size())
	{
		delete[]Matrix::pVector;
		Matrix::size_ = in_list.size();
		Matrix::pVector = new Vector<ValType>[Matrix::size_];
	}

	auto it = in_list.begin();
	for (int i = 0; i < Matrix::size_; i++)
		Matrix::pVector[i] = *it, it++;

	return *this;
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator= (const Vector<ValType>& Vec)
{
	for (int i = 0; i < this->size_; i++)
		this->pVector[i][col_ind] = Vec[i];
	return *this;
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator= (std::initializer_list<ValType> init_list)
{
	auto it = init_list.begin();
	for (int i = 0; i < this->size_; i++)
		this->pVector[i][col_ind] = *it, it++;
	return *this;
}

/*

	Операции с матрицами

*/

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator*(const Matrix<ValType>& mt)
{
	int row = Matrix::size_, cell = mt.pVector[0].size(), n = mt.size_;
	Matrix<ValType> res(row, cell);
	for (int i = 0; i < row; i++)
		for (int j = 0; j < cell; j++)
			for (int k = 0; k < n; k++)
				res.pVector[i][j] += Matrix::pVector[i][k] * mt.pVector[k][j];
	return res;
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator+(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator+(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator+=(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator+=(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator-=(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator-=(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator-(const Matrix<ValType>& mt)
{
	return Vector<Vector<ValType>>::operator-(mt);
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator*=(const Matrix<ValType>& mt)
{
	Matrix<ValType> res = *this;
	*this = *this * mt;
	return res;
}

/*

	Полезные операции

*/

template <class ValType>
Matrix<ValType> Matrix<ValType>::t()
{
	return TransposedMatrix(*this);
}

template <class ValType>
Matrix<ValType> Matrix<ValType>::i()
{
	return InverseMatrix(*this);
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::diag(const Vector<ValType>& diag_el)
{
	if (this->size_ != diag_el.size_)
		*this = Matrix<ValType>(diag_el.size_);
	for (int i = 0; i < this->size_; i++)
		this->pVector[i][i] = diag_el[i];
	return *this;
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::diag(std::initializer_list<ValType> init_list)
{
	if (this->size_ != init_list.size())
		*this = Matrix<ValType>(init_list.size());

	auto it = init_list.begin();
	for (int i = 0; i < this->size_; i++)
		this->pVector[i][i] = *it, it++;
	return *this;
}

/*

	Операции с числами

*/

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator*(const ValType& Val)
{
	Matrix<ValType> res = *this;
	for (int i = 0; i < res.size_; i++)
		res[i] = Matrix::pVector[i] * Val;
	return res;
}

template <class ValType>
Matrix<ValType>  Matrix<ValType>::operator/(const ValType& Val)
{
	Matrix<ValType> res = *this;
	for (int i = 0; i < res.size_; i++)
		res[i] = Matrix::pVector[i] / Val;
	return res;
}

/*

	Логические операции

*/

template <class ValType>
bool Matrix<ValType>::operator==(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator==(mt);
}

template <class ValType>
bool Matrix<ValType>::operator!=(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator!=(mt);
}

template <class ValType>
bool Matrix<ValType>::operator>(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator>(mt);
}

template <class ValType>
bool Matrix<ValType>::operator<(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator<(mt);
}

template <class ValType>
bool Matrix<ValType>::operator>=(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator>=(mt);
}

template <class ValType>
bool Matrix<ValType>::operator<=(const Matrix<ValType>& mt) const
{
	return Vector<Vector<ValType>>::operator<=(mt);
}

/*

	Инкремент и декремент

*/

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator++()
{
	*this = Vector<Vector<ValType>>::operator++();
	return *this;
}

template <class ValType>
Matrix<ValType> Matrix<ValType>::operator++(int)
{
	return Vector<Vector<ValType>>::operator++(0);
}

template <class ValType>
Matrix<ValType>& Matrix<ValType>::operator--()
{
	*this = Vector<Vector<ValType>>::operator--();
	return *this;
}

template <class ValType>
Matrix<ValType> Matrix<ValType>::operator--(int)
{
	return Vector<Vector<ValType>>::operator--(0);
}

/*

	Полезно

*/

template <class ValType>
Vector<ValType> operator*(const Vector<ValType>& Vec, const Matrix<ValType>& mt)
{
	int n = Vec.size(), cell = mt[0].size();
	Vector<ValType> res(cell);
	for (int i = 0; i < cell; i++)
		for (int j = 0; j < n; j++)
			res[i] += Vec[j] * mt[j][i];
	return res;
}

template <class ValType>
Vector<ValType> operator*(const Matrix<ValType>& mt, const Vector<ValType>& Vec)
{
	Vector<ValType> res(mt.size());
	for (int i = 0; i < mt.size(); i++)
		res[i] = mt[i] * Vec;
	return res;
}

template <class ValType>
Matrix<ValType> Single(int size_)
{
	Matrix<ValType> res(size_);
	for (int i = 0; i < res.size(); i++)
		res[i][i] = 1;
	return res;
}

template <class ValType>
Matrix<ValType> pow(const Matrix<ValType>& mt, int degree)
{
	Matrix<ValType> res = Single<ValType>(mt.size());
	for (int i = 1; i <= degree; i++)
		res *= mt;
	return res;
}

template <class ValType>
Matrix<ValType> InverseMatrix(const Matrix<ValType>& mt)
{
	Matrix<ValType> res = Single<ValType>(mt.size()), matr = mt;
	ValType sup;
	for (int i = 1; i < matr.size(); i++)
	{
		if (fabs(matr[i - 1][i - 1]) <= pow(10, -14))
			for (int j = i; j < matr.size(); j++)
			{
				matr[i - 1] = matr[i - 1] + matr[j];
				res[i - 1] = res[i - 1] + res[j];
			}
		for (int j = i; j < matr.size(); j++)
		{
			sup = matr[j][i - 1] / matr[i - 1][i - 1];
			matr[j] = matr[j] - matr[i - 1] * sup;
			res[j] = res[j] - res[i - 1] * sup;
		}
	}
	for (int i = matr.size() - 2; i >= 0; i--)
	{
		for (int j = i; j >= 0; j--)
		{
			sup = matr[j][i + 1] / matr[i + 1][i + 1];
			matr[j] = matr[j] - matr[i + 1] * sup;
			res[j] = res[j] - res[i + 1] * sup;
		}
	}
	for (int i = 0; i < matr.size(); i++)
		res[i] = res[i] / matr[i][i];
	return res;
}

template <class ValType>
Matrix<ValType> TransposedMatrix(const Matrix<ValType>& mt)
{
	int row = mt.size(), cell = mt[0].size();
	Matrix<ValType> res(cell, row);
	for (int i = 0; i < cell; i++)
		for (int j = 0; j < row; j++)
			res[i][j] = mt[j][i];
	return res;
}

template <class ValType>
ValType Determinant(const Matrix<ValType>& mt)
{
	ValType res = ValType();

	if (mt.size() == 1)
		return mt[0][0];
	if (mt.size() == 2)
		return mt[0][0] * mt[1][1] - mt[0][1] * mt[1][0];

	Matrix<ValType> SupMT(mt.size() - 1);
	for (int i = 0; i < mt.size(); i++)
	{
		SupMT = Minor(mt, 0, i);
		res += pow(-1, i) * Determinant(SupMT) * mt[0][i];
	}

	return res;
}

template <class ValType>
Matrix<ValType> Minor(const Matrix<ValType>& mt, int ntrow, int ntcell)
{
	int row = mt.size(), cell = mt[0].size();
	Matrix<ValType> res(row - 1, cell - 1);
	for (int i = 0, ires = 0; i < row; i++)
	{
		if (i != ntrow)
		{
			for (int j = 0, jres = 0; j < cell; j++)
				if (j != ntcell)
					res[ires][jres] = mt[i][j], jres++;
			ires++;
		}
	}
	return res;
}

template <class ValType>
Matrix<ValType> CellRowMultiply(const Vector<ValType>& cell, const Vector<ValType>& row)
{
	Matrix<ValType> res(cell.size(), row.size());
	for (int i = 0; i < cell.size(); i++)
		for (int j = 0; j < row.size(); j++)
			res[i][j] = cell[i] * row[j];
	return res;
}

template <class ValType>
Matrix<ValType> operator*(const ValType& Val, const Matrix<ValType>& mt)
{
	Matrix<ValType> res(mt.size());
	for (int i = 0; i < mt.size(); i++)
		res[i] = Val * mt[i];
	return res;
}

#endif