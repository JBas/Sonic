#include <stdint.h>
#include <ncurses.h>
// #include <linux/spi/spidev.h>

// #define SPI_DEVICE "/dev/spidev0.0"

struct update_header {
    uint16_t x;
    uint16_t y;
    uint16_t z;
    uint8_t r;
    uint8_t g;
    uint8_t b;
} typedef update_t;

struct frame_header {
    update_t *updates;
} typedef frame_t;

void send_update(update_t* u) {
    return;
}

void reset() {
    update_t u = { 0 }; // (0, 0, 0) represents the middle of the display
    send_update(&u);
}


int main() {

    initscr();

    geti

    endwin();

    reset();

    return 0;

}