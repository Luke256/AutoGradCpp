# AutoGradCpp

C++用簡易自動微分ヘッダー(リバースモード)

# 使用方法

`AutoGrad.hpp`を適当な場所に配置し、プログラム内でインクルードすれば使用できます。

# 注意

このヘッダーは簡易的な自動微分を実装しています。つまり、**最後の計算結果についてのみ**正しく微分が行われることを保証します。計算の途中経過に対して`backward()`メソッドを使用しても、正しい微分が得られるかどうかは未検証です。

# リファレンス
## AutoGrad::Varクラス
変数をVarクラスで定義すると自動微分を行えるようになります。

例：
```cpp
AutoGrad::Var x(100);
```

## AutoGrad::Var::value() -> double&
格納されている値の情報を返します

## AutoGrad::Var::grad() -> double&
格納されている勾配の情報を返します。微分がされていない場合、0を返します。

## AutoGrad::Var::backward() -> void
偏微分を実行します。

> [!Warning]
> この関数を呼び出す前に、計算に使用したすべての変数の勾配を0で初期化しなければならないことに注意してください。

例：
```cpp
AutoGrad::Var a(1), b(2);
AutoGrad::Var c = a * b;

c.backward();

a.grad() // 2
b.grad() // 1
```

## AutoGrad::exp(const Var& v) -> Var
指数関数を適用させた値を返します。

## AutoGrad::log(const Var& v) -> Var
自然対数を適用させた値を返します。

## AutoGrad::sin(const Var& v) -> Var
正弦関数を適用させた値を返します。

## AutoGrad::cos(const Var& v) -> Var
余弦関数を適用させた値を返します。

## AutoGrad::tan(const Var& v) -> Var
正接関数を適用させた値を返します。

## AutoGrad::ClearGradientTape() -> void
計算の履歴を消去します。

計算履歴が不要になった時に呼び出してください。