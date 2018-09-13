using namespace std;

class Shape {
    protected:
        float ** points_;
        unsigned int npoints_;

    public:
        Shape();
        ~Shape();
        virtual float** get() const = 0;
};