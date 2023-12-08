#include "cell.cpp"
#include "random.cpp"
#include <vector>
#include <fstream>
#include <cmath>
#include <math.h>
#include <bits/stdc++.h> //for sorting
template <class T>
/*
DONE:
- [x] Szabo
- [x] addition of cells Conditional  on presence of others
- [x] Connect the boundary cells
*/

/*
TODO:
- [x] Make the source an arc (instead of a line)
- [x] Add repulsion from arc
- [] Find the cell closest to the arc than the top one
- [] Set distance threshold on boundarySequence construction
- [] Add attraction between boundary cells
- [] Fix long distance attraction
*/
class Simulation
{
private:
public:
    std::vector<Cell<T>> cells;
    std::vector<Cell<T> *> neighborCells;
    std::vector<Cell<T> *> boundaryCells;
    std::vector<Cell<T> *> bNeighborCells;
    std::vector<Cell<T> *> validBoundaryCells;
    std::vector<Cell<T> *> boundaryCellsSequence;
    std::ofstream outfile;
    std::ofstream boundaryOutfile;
    std::ofstream validBoundaryOutfile;
    std::ofstream boundarySequenceOutfile;
    T deltat;
    T enterGapT, enterT;
    T width, height, G, gmax;
    T rt;
    T brt, fov; // boundary detection thresh
    T tmax;
    T beta;
    T FMaxIn;
    T progt, progtstep, progtmax, nprogsteps;
    T alpha;
    T ythresh;
    VEC2<T> *circleCenter;
    T maxCircleRadius;

    int N;

    Simulation(int N = 70,
               T tmax = 100, T deltat = 0.01,
               T width = 1200, T height = 1200)
    {
        // beta = 60;
        brt = 300;
        fov = M_PI_2;
        enterT = 0;
        enterGapT = deltat;
        beta = 30;
        FMaxIn = 50;
        alpha = 0.005;
        ythresh = height - 50;
        // G = 30;
        circleCenter = new VEC2<T>(width / 2, height + 100);
        maxCircleRadius = 185;
        // gmax = 10;
        rt = 70;
        outfile.open("outdata.csv");
        boundaryOutfile.open("boundaryOutData.csv");
        validBoundaryOutfile.open("topBoundaryOutData.csv");
        boundarySequenceOutfile.open("boundarySequenceOutData.csv");
        this->deltat = deltat;
        this->tmax = tmax;
        this->N = N;
        this->width = width;
        this->height = height;
        setupProgressBar();

        // for (int i = 0; i < N; i++)
        // {
        //     Cell<T> cell;
        //     cell.setRandomPosition(width, height);
        //     cells.push_back(cell);
        //     VEC2<T> f;
        //     forces.push_back(f);
        // }
    }
    void addCell(T t)
    {
        if (t >= enterT)
        {
            enterT += enterGapT;
            Cell<T> cell;
            cell.pos.x = Rand<T>::getUniformRand(width / 3.5, width - width / 3.5);
            cell.vel.y = Rand<T>::getUniformRand(2, 4);
            cell.vel.x = Rand<T>::getUniformRand(-0.5, 0.5);
            getNeighbors(&cell, rt * 0.7);
            if (neighborCells.size() == 0)
                cells.push_back(cell);
            else
                cell.~Cell();
        }
    }
    void addCellArc(T t, VEC2<T> *circleCenter, T maxRadius)
    {
        if (t >= enterT)
        {
            int nTries = 10;
            while (nTries-- > 0)
            {
                enterT += enterGapT;
                Cell<T> cell;
                T theta = Rand<T>::getUniformRand(0, M_PI * 2);
                T length = Rand<T>::getUniformRand(maxRadius - 60, maxRadius);
                VEC2<T> unit(cos(theta), sin(theta));
                VEC2<T> pos = *circleCenter + unit * length;
                if (pos.y > height)
                    return;
                T vmag = Rand<T>::getUniformRand(1, 3);
                VEC2<T> vel = unit * vmag;
                // cell.vel.y = Rand<T>::getUniformRand(2, 4) * (-1);
                cell.pos = pos;
                cell.vel = vel;
                getNeighbors(&cell, rt * 0.4);
                if (neighborCells.size() == 0)
                    cells.push_back(cell);
                else
                    cell.~Cell();
            }
        }
    }
    void simulate()
    {
        for (T t = 0; t < tmax; t += deltat)
        {
            // addCell(t);
            addCellArc(t, circleCenter, maxCircleRadius);
            update();
            printPosToFile();
            printBoundaryPosToFile();
            printValidBoundaryPosToFile();
            printBoundarySequenceToFile();
            updateProgressBar(t);
        }
        std::cout << "Done" << std::endl;
    }
    void update()
    {
        // for (int i = 0; i < N; i++)
        for (auto &cell : cells)
        {
            getNeighbors(&cell, rt);
            VEC2<T> FIn = getInteractionForceSzabo(&cell);
            VEC2<T> FCircRep = getCircleRepulsiveForce(&cell);
            // VEC2<T> FIn = getInteractionForce(&cell);
            VEC2<T> FVc = getVicsekForce(&cell);
            VEC2<T> FVisc = cell.vel * (-alpha);

            // VEC2<T> FNo = getNoiseForce(&cell);
            // if (FIn.mag() > FMaxIn){
            //     FIn = FIn / FIn.mag();
            //     FIn = FIn * FMaxIn;
            // }

            cell.acc += FCircRep;
            // cell.acc += FIn;
            cell.pos += FIn * deltat; // only for Szabo
            cell.pos += VEC2<T>(Rand<T>::getUniformRand(-1, 1),
                                Rand<T>::getUniformRand(-1, 1)) *
                        0.3;
            // cell.acc += FVc * beta;
            // cell.acc += FNo;
            cell.acc += FVisc;
        }
        calcBoundaryCells(brt, fov);
        calcBoundaryCellsOutsideCircle();
        calcBoundaryCellsSequence(validBoundaryCells);
        if (boundaryCellsSequence.size() >= 3)
        {
            for (int i = 1; i < boundaryCellsSequence.size() - 1; i++)
            {
                Cell<T> *current = boundaryCellsSequence.at(i);
                Cell<T> *previous = boundaryCellsSequence.at(i - 1);
                Cell<T> *next = boundaryCellsSequence.at(i + 1);

                VEC2<T> FBndIn = getInteractionForceSzabo(current, previous) +
                                 getInteractionForceSzabo(current, next);
                FBndIn = FBndIn * 2;
                current->pos += FBndIn * deltat;
            }
        }
        for (auto &cell : cells)
        {
            cell.update(deltat);
        }
    }
    VEC2<T> getInteractionForceSzabo(Cell<T> *A, Cell<T> *B)
    {
        VEC2<T> D(B->pos - A->pos);
        T r = D.mag();
        T ifMag = interactionForceSzabo(r);
        return D.unit() * ifMag;
    }
    VEC2<T> getInteractionForceSzabo(Cell<T> *cell)
    {
        VEC2<T> F;
        for (auto neighborCell : neighborCells)
        {
            F += getInteractionForceSzabo(cell, neighborCell);
        }
        return F;
    }
    // VEC2<T> getInteractionForce(Cell<T> *cell)
    // {
    //     VEC2<T> F;
    //     for (auto neighborCell : neighborCells)
    //     {
    //         Cell<T> A = *cell;
    //         Cell<T> B = *neighborCell;
    //         VEC2<T> D(B.pos - A.pos);
    //         T r = D.mag();
    //         // T ifMag = interactionForceSpringLike(r);
    //         // T ifMag = interactionForceSzabo(r);
    //         T ifMag = interactionForceMagNirgov(r);
    //         VEC2<T> iForce = D.unit() * ifMag;
    //         F += iForce;
    //     }
    //     return F;
    // }
    VEC2<T> getCircleRepulsiveForce(Cell<T> *cell)
    {
        VEC2<T> D = cell->pos - *circleCenter;
        T r = D.mag();
        if (r > maxCircleRadius)
            return VEC2<T>(0, 0);

        T repmag = 100 * (1 / r);
        return D.unit() * repmag;
    }
    VEC2<T> getVicsekForce(Cell<T> *cell)
    {
        VEC2<T> FVc;
        int ninr = 0;
        for (auto neighborCell : neighborCells)
        {
            VEC2<T> dv = neighborCell->vel - cell->vel;
            FVc += dv;
            ninr++;
        }
        if (ninr == 1)
        { // no neighbors except self
            return VEC2<T>(0, 0);
        }
        else
        {
            return FVc / T(ninr);
        }
    }
    VEC2<T> unitnoise()
    {
        VEC2<T> U = VEC2<T>(Rand<T>::getUnitRand(), Rand<T>::getUnitRand());
        VEC2<T> xi(sqrt(-2 * log(U.x)) * cos(2 * M_PI * U.y),
                   sqrt(-2 * log(U.x)) * sin(2 * M_PI * U.y));
        return xi;
    }
    T getLocalCellDensity()
    {
        T rho = neighborCells.size() + 1;
        return rho / (M_PI * rt * rt);
    }
    VEC2<T> getNoiseForce(Cell<T> *cell)
    {
        T sig0 = 150.,
          sig1 = 300.,
          rho0 = N / (width * height),
          rho = 0.,
          sig = 0.,
          tau = 1.39;
        VEC2<T> xi = unitnoise();
        rho = getLocalCellDensity();
        sig = sig0 + (sig1 - sig0) * (1 - rho / rho0);
        // integrate and normalize
        cell->eta -= cell->eta * deltat / tau;
        cell->eta += xi * sqrt(deltat) / tau;
        cell->eta = cell->eta / cell->eta.mag();
        return cell->eta * sig;
    }
    void printPosToFile()
    {
        for (auto cell : cells)
        {
            outfile << cell.pos.x
                    << " "
                    << cell.pos.y
                    << " ";
        }
        outfile << std::endl;
    }
    void printBoundaryPosToFile()
    {
        for (auto cell : boundaryCells)
        {

            boundaryOutfile << cell->pos.x
                            << " "
                            << cell->pos.y
                            << " ";
        }
        boundaryOutfile << std::endl;
    }
    void printValidBoundaryPosToFile()
    {
        for (auto cell : validBoundaryCells)
        {

            validBoundaryOutfile << cell->pos.x
                                 << " "
                                 << cell->pos.y
                                 << " ";
        }
        validBoundaryOutfile << std::endl;
    }
    void printBoundarySequenceToFile()
    {
        for (auto cell : boundaryCellsSequence)
        {

            boundarySequenceOutfile << cell->pos.x
                                    << " "
                                    << cell->pos.y
                                    << " ";
        }
        boundarySequenceOutfile << std::endl;
    }
    void getNeighbors(Cell<T> *cell, T rt)
    /*Get list of neighbors. EXLUDES self*/
    {
        // reset neighborCells
        neighborCells.clear();
        for (int i = 0; i < cells.size(); i++)
        {
            if (cell != &cells.at(i))
            {
                VEC2<T> D = cells.at(i).pos - cell->pos;
                T r = D.mag();
                if (r < rt)
                {
                    neighborCells.push_back(&cells.at(i));
                }
            }
        }
    }
    void getbNeighbors(Cell<T> *cell, std::vector<Cell<T> *> &bcells, T rt)
    {
        // reset neighborCells
        bNeighborCells.clear();
        for (int i = 0; i < bcells.size(); i++)
        {
            if (cell != bcells.at(i))
            {
                VEC2<T> D = bcells.at(i)->pos - cell->pos;
                T r = D.mag();
                if (r < rt)
                {
                    bNeighborCells.push_back(bcells.at(i));
                }
            }
        }
    }
    void calcBoundaryCells(T brt, T fov)
    {
        boundaryCells.clear();
        for (auto &cell : cells)
        {
            getNeighbors(&cell, brt);
            if (isOnBoundary(&cell, brt, fov))
                boundaryCells.push_back(&cell);
        }
    }
    void calcBoundaryCellsBelowY(T ythresh)
    { // below y-threshold
        validBoundaryCells.clear();
        for (auto &cell : boundaryCells)
        {
            if (cell->pos.y < ythresh)
                validBoundaryCells.push_back(cell);
        }
    }
    void calcBoundaryCellsOutsideCircle()
    { // outside circle
        validBoundaryCells.clear();
        for (auto &cell : boundaryCells)
        {
            if ((cell->pos - *circleCenter).mag() > maxCircleRadius)
                validBoundaryCells.push_back(cell);
        }
    }
    Cell<T> *gettopCell(std::vector<Cell<T> *> &bcells)
    {
        Cell<T> *topCell;
        topCell = bcells.at(0);
        for (auto &cell : bcells)
        {
            if (cell->pos.x > topCell->pos.x)
            {
                topCell = cell;
            }
        }
        return topCell;
    }
    Cell<T> *getNearestCell(Cell<T> *currentCell, std::vector<Cell<T> *> bNeighborCells)
    {
        T shortest = FLT_MAX;
        Cell<T> *nearestCell;
        for (auto &cell : bNeighborCells)
        {
            if (currentCell != cell) // redundant
            {
                T dist = (currentCell->pos - cell->pos).mag();
                if (dist < shortest)
                {
                    shortest = dist;
                    nearestCell = cell;
                }
            }
        }
        return nearestCell;
    }
    // void calcBoundaryCellsSequence(std::vector<Cell<T> *> bcells)
    // {
    //     boundaryCellsSequence.clear();
    //     if (bcells.size() >= 10)
    //     {
    //         Cell<T> *topCell = getBottomMostCell(bcells);
    //         boundaryCellsSequence.push_back(topCell);
    //         Cell<T> *currentCell = topCell;
    //         while (bcells.size() > 2) //why should it be 2???
    //         {
    //             getbNeighbors(currentCell, bcells, brt);
    //             Cell<T> *nearestCell = getNearestCell(currentCell);
    //             boundaryCellsSequence.push_back(nearestCell);
    //             auto it = std::find(bcells.begin(), bcells.end(), currentCell);
    //             currentCell = nearestCell;
    //             if (it != bcells.end())
    //             {
    //                 bcells.erase(it);
    //             }
    //         }
    //     }
    // }
    Cell<T> *getBottomMostCell(std::vector<Cell<T> *> cells)
    {
        Cell<T> *bottomMostCell = cells.at(0);
        for (auto &cell : cells)
        {
            if (cell->pos.y < bottomMostCell->pos.y)
            {
                bottomMostCell = cell;
            }
        }
        return bottomMostCell;
    }
    // void calcBoundaryCellsSequence(std::vector<Cell<T> *> bcells)
    // {
    //     boundaryCellsSequence.clear();
    //     if (bcells.size() > 10)
    //     {
    //         Cell<T> *bottomMostCell = getBottomMostCell(bcells);
    //         boundaryCellsSequence.push_back(bottomMostCell);
    //         Cell<T> *currentCell = bottomMostCell;
    //         // iterate
    //         while (bcells.size() > 1)
    //         {
    //             // find nearest cell to current cell in bcells
    //             getbNeighbors(currentCell, bcells, brt);
    //             Cell<T> *nearestCell = getNearestCell(currentCell);
    //             boundaryCellsSequence.push_back(nearestCell);
    //             removeElement(bcells, currentCell);
    //             currentCell = nearestCell;
    //             // add nearest cell
    //             // remove current cell
    //             // set nearest cell as current cell
    //         }
    //     }
    // }

    // void calcBoundaryCellsSequence(std::vector<Cell<T> *> bcells)
    // { // semi-final working
    //     boundaryCellsSequence.clear();
    //     if (bcells.size() > 10)
    //     {
    //         Cell<T> *topCell = getBottomMostCell(bcells);
    //         boundaryCellsSequence.push_back(topCell);
    //         Cell<T> *currentCell = topCell;
    //         removeElement(bcells, currentCell);
    //         while (bcells.size() > 0)
    //         {
    //             getbNeighbors(currentCell, bcells, width);
    //             Cell<T> *nearestCell = getNearestCell(currentCell);
    //             boundaryCellsSequence.push_back(nearestCell);
    //             currentCell = nearestCell;
    //             removeElement(bcells, currentCell);
    //         }
    //     }
    // }
    void calcBoundaryCellsSequence(std::vector<Cell<T> *> bcells)
    { // using new functions
        boundaryCellsSequence.clear();
        if (bcells.size() > 10)
        { // find top most cell
            auto topCell = gettopCell(bcells);
            boundaryCellsSequence.push_back(topCell);
            removeElement(bcells, topCell);
            auto currentCell = topCell;
            // find closest cell
            while (bcells.size() > 0)
            {
                getbNeighbors(currentCell, bcells, 600);
                if (bNeighborCells.size() == 0)
                    break;
                auto nearestCell = getNearestCell(currentCell, bNeighborCells);
                boundaryCellsSequence.push_back(nearestCell);
                removeElement(bcells, nearestCell);
                currentCell = nearestCell;
            }
        }
    }
    void calcBoundaryCellsSequencea(std::vector<Cell<T> *> bcells)
    { // manual, without previous functions
        boundaryCellsSequence.clear();
        if (bcells.size() > 10)
        { // find top most cell
            auto topCell = bcells.at(0);
            for (auto &cell : bcells)
            {
                if (cell->pos.y > topCell->pos.y)
                    topCell = cell;
            }
            boundaryCellsSequence.push_back(topCell);
            removeElement(bcells, topCell);
            auto currentCell = topCell;
            // find closest cell
            while (bcells.size() > 0)
            {
                auto nearestCell = bcells.at(0);
                T shortestDist = FLT_MAX;
                for (auto &bcell : bcells)
                {
                    if (currentCell != bcell)
                    {
                        T dist = (currentCell->pos - bcell->pos).mag();
                        if (dist < shortestDist)
                        {
                            shortestDist = dist;
                            nearestCell = bcell;
                        }
                    }
                }
                boundaryCellsSequence.push_back(nearestCell);
                removeElement(bcells, nearestCell);
                currentCell = nearestCell;
            }
        }
    }

    void removeElement(std::vector<Cell<T> *> &bcells, Cell<T> *currentCell)
    {
        for (int i = 0; i < bcells.size(); i++)
        {
            if (bcells.at(i) == currentCell)
            {
                bcells.erase(bcells.begin() + i);
                // return;
            }
        }
    }
    bool isOnBoundary(Cell<T> *cell, T brt, T fov)
    {
        int nelt = neighborCells.size();
        if (nelt <= 1)
            return true;
        std::vector<T> q;
        for (auto &neighborCell : neighborCells)
        {
            VEC2<T> s = neighborCell->pos - cell->pos;
            q.push_back(atan2(s.y, s.x));
        }
        std::sort(q.begin(), q.end()); // ascending
        T lg = findLargestGap(q);
        return (lg > fov);
    }
    T findLargestGap(const std::vector<T> &q)
    {
        T lg = 0.;
        T t;
        int nelt = q.size();
        for (int i = 0; i < nelt - 1; i++)
        {
            t = q.at(i + 1) - q.at(i);
            if (t > lg)
            {
                lg = t;
            }
        }
        t = 2 * M_PI - q.at(nelt - 1) + q.at(0);
        if (t > lg)
        {
            lg = t;
        }
        return lg;
    }
    T interactionForceMagNirgov(T r) // nirgov
    {
        if (r < rt)
        {
            T U0 = 2650, U1 = 30, U2 = 2, U3 = 1;
            T A0 = 8, A1 = 2, A2 = 25, A3 = 26;
            T force = 0;
            force += U0 * r * exp(-(pow((r / A0), 2)));
            force += U2 * exp(-r / A2);
            force -= U3 * pow(r - A3, 2) * Hv(r - A3);
            force += U1 * (r - A1) * Hv(r - A1);
            return -force;
        }
        else
        {
            return 0.;
        }
    }
    // T interactionForceSpringLike(T r)
    // {
    //     T fmag = 1 / pow(r - 5, 2) + 0.15 * r - 10.0;
    //     return fmag;
    // }
    T interactionForceSzabo(T r)
    {

        T Fadh = -10.0;
        T Frep = 10.0;

        T Req = 60.0;
        T R0 = 100.0;
        T Fmag = 0.0;
        if (r < Req)
            Fmag = Frep * (Req - r) / Req;
        else
            Fmag = Fadh * (r - Req) / (R0 - Req);

        return -Fmag;
    }
    T interactionForceSzaboBoundary(T r)
    {

        T Fadh = -50.0;
        T Frep = 10.0;

        T Req = 60.0;
        T R0 = 100.0;
        T Fmag = 0.0;
        if (r < Req)
            Fmag = Frep * (Req - r) / Req;
        else
            Fmag = Fadh * (r - Req) / (R0 - Req);

        return -Fmag;
    }

    T gravityForceMag(T r)
    {
        return (G / r);
    }
    T Hv(T r)
    {
        return T(r > 0);
    }
    void setupProgressBar()
    {
        progt = 0;
        progtmax = tmax;
        nprogsteps = 30;
        progtstep = progtmax / nprogsteps;
        std::cout << "|";
        for (int i = 0; i < nprogsteps; i++)
        {
            std::cout << "-";
        }
        std::cout << "|" << std::endl;
        std::cout << "|=";
    }
    void updateProgressBar(T t)
    {
        if (t > progt)
        {
            progt += progtstep;
            if (progt >= tmax - progtstep + 0.01)
                std::cout << "|" << std::endl;
            else
                std::cout << "=" << std::flush;
        }
    }
};
