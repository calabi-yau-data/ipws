#ifndef SETTINGS_H
#define SETTINGS_H

struct Settings {
    unsigned redundancy_check_skip_recursions{};
    bool allow_weights_one{};
    bool allow_weights_one_half{};
    bool allow_weights_sum_one{};
    bool print_statistics{};
    bool print_weight_systems{};
    bool debug_ignore_symmetries{};
    bool debug_disable_lex_order{};
};

extern Settings g_settings;

#endif
