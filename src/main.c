#include <constraint.h>
#include <globaldef.h>
#include <simulation.h>
#include <unistd.h>

int main() {
    simulation *sim = init_simulation(3, 3, 5);
    if (!sim)
        return 1;

    for (int i = 0; i < 10; i++) {
        propagate_simulation(sim, .01);
        output_positions(sim);
        usleep(1e6);
    }
    destruct_simulation(sim);
    return 0;
}