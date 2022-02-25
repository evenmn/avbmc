#include <iostream>
#include <string>
#include <vector>


class Constraint
{
public:
    Constraint(class Box *);
    virtual bool verify() = 0;

    std::string label;

protected:
    class Box* box = nullptr;
};
