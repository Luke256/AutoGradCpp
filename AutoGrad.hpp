# include <iostream>
# include <memory>
# include <vector>
# include <functional>
# include <cmath>

namespace AutoGrad
{

using _BFty = std::function<void(std::shared_ptr<double>,std::shared_ptr<double>,std::shared_ptr<double>,std::shared_ptr<double>,std::shared_ptr<double>)>;

namespace Backwards
{

    void Add(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg;
        *bg += *cg;
    }

    void Sub(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg;
        *bg += - *cg;
    }

    void Mul(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg * *bv;
        *bg += *cg * *av;
    }

    void Div(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg / *bv;
        *bg += - *cg * *av / (*bv * *bv);
    }

    void Exp(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg * std::exp(*av);
    }

    void Log(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg / *av;
    }

    void Sin(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg * std::cos(*av);
    }

    void Cos(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg * -std::sin(*av);
    }

    void Tan(std::shared_ptr<double>av, std::shared_ptr<double>bv, std::shared_ptr<double>ag, std::shared_ptr<double>bg, std::shared_ptr<double>cg)
    {
        *ag += *cg * (1 + std::pow(std::tan(*av), 2));
    }
}

struct History
{
    std::shared_ptr<double>av,bv,ag,bg,cg;
    _BFty backwardFunction;
    void backward()
    {
        backwardFunction(av, bv, ag, bg, cg);
    }
};

struct Var
{
    std::shared_ptr<double> value_ptr, grad_ptr;
    int tapePosition;

    Var(): 
        value_ptr(std::make_shared<double>(0)), 
        grad_ptr(std::make_shared<double>(0)) 
    {

    }

    Var(double v) : 
        value_ptr(std::make_shared<double>(v)), 
        grad_ptr(std::make_shared<double>(0)) 
    {

    }
    
    static std::vector<History>& GRADIENT_TAPE()
    {
        static std::vector<History> tape;
        return tape;
    }

    double& value()
    {
        return *value_ptr;
    }

    double& grad()
    {
        return *grad_ptr;
    }

    void backward()
    {
        *grad_ptr = 1;
        for (int i = tapePosition; i >= 0; i--)
        {
            GRADIENT_TAPE()[i].backward();
        }
    }

    Var& operator+=(const Var& v)
    {
        Var next;
        *next.value_ptr = *value_ptr + *v.value_ptr;
        next.tapePosition = GRADIENT_TAPE().size();
        GRADIENT_TAPE().push_back({value_ptr, v.value_ptr, grad_ptr, v.grad_ptr, next.grad_ptr, Backwards::Add});
        *this = next;
        return *this;
    }

    Var& operator-=(const Var& v)
    {
        Var next;
        *next.value_ptr = *value_ptr - *v.value_ptr;
        next.tapePosition = GRADIENT_TAPE().size();
        GRADIENT_TAPE().push_back({value_ptr, v.value_ptr, grad_ptr, v.grad_ptr, next.grad_ptr, Backwards::Sub});
        *this = next;
        return *this;
    }

    Var& operator*=(const Var& v)
    {
        Var next;
        *next.value_ptr = *value_ptr * *v.value_ptr;
        next.tapePosition = GRADIENT_TAPE().size();
        GRADIENT_TAPE().push_back({value_ptr, v.value_ptr, grad_ptr, v.grad_ptr, next.grad_ptr, Backwards::Mul});
        *this = next;
        return *this;
    }

    Var& operator/=(const Var& v)
    {
        Var next;
        *next.value_ptr = *value_ptr / *v.value_ptr;
        next.tapePosition = GRADIENT_TAPE().size();
        GRADIENT_TAPE().push_back({value_ptr, v.value_ptr, grad_ptr, v.grad_ptr, next.grad_ptr, Backwards::Div});
        *this = next;
        return *this;
    }

    Var operator+(const Var& v) { return Var{*this} += v; }
    Var operator-(const Var& v) { return Var{*this} -= v; }
    Var operator*(const Var& v) { return Var{*this} *= v; }
    Var operator/(const Var& v) { return Var{*this} /= v; }
};

void ClearGradientTape()
{
    Var::GRADIENT_TAPE().clear();
}

Var exp(Var& v)
{
    Var res{ *v.value_ptr };
    *res.value_ptr = std::exp(*res.value_ptr);
    res.tapePosition = Var::GRADIENT_TAPE().size();
    Var::GRADIENT_TAPE().push_back({v.value_ptr, 0, v.grad_ptr, 0, res.grad_ptr, Backwards::Exp});
    return res;
}

Var log(Var& v)
{
    Var res{ *v.value_ptr };
    *res.value_ptr = std::log(*res.value_ptr);
    res.tapePosition = Var::GRADIENT_TAPE().size();
    Var::GRADIENT_TAPE().push_back({v.value_ptr, 0, v.grad_ptr, 0, res.grad_ptr, Backwards::Log});
    return res;
}

Var sin(Var& v)
{
    Var res{ *v.value_ptr };
    *res.value_ptr = std::sin(*res.value_ptr);
    res.tapePosition = Var::GRADIENT_TAPE().size();
    Var::GRADIENT_TAPE().push_back({ v.value_ptr, 0, v.grad_ptr, 0, res.grad_ptr, Backwards::Sin });
}

Var cos(Var& v)
{
    Var res{ *v.value_ptr };
    *res.value_ptr = std::cos(*res.value_ptr);
    res.tapePosition = Var::GRADIENT_TAPE().size();
    Var::GRADIENT_TAPE().push_back({ v.value_ptr, 0, v.grad_ptr, 0, res.grad_ptr, Backwards::Cos });
}

Var tan(Var& v)
{
    Var res{ *v.value_ptr };
    *res.value_ptr = std::tan(*res.value_ptr);
    res.tapePosition = Var::GRADIENT_TAPE().size();
    Var::GRADIENT_TAPE().push_back({ v.value_ptr, 0, v.grad_ptr, 0, res.grad_ptr, Backwards::Tan });
}


}