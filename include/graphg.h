#ifndef GRAPHG_H
#define GRAPHG_H
#include <vector>

using namespace std;


class Node
{
    public:
        Node(int _id);
        virtual ~Node();
        int id;
        float f; //fixed cost
        float p; // profit

    protected:

    private:
};

#endif // GRAPHG_H
