#include <gtk/gtk.h>


void
on_amin_toggled                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
expose_gl                              (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_amax_toggled                        (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_spinbuttonmin_changed               (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_spinbuttonmax_changed               (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_open1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_save1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_export1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_scale_changed                       (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_delete_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_color_clicked                       (GtkButton       *button,
                                        gpointer         user_data);

void
on_constant1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_flat1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_smooth1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_csmooth1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_color_selector_clicked              (GtkButton       *button,
                                        gpointer         user_data);

void
on_properties1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_close1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_clear1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_properties1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_preferences1_activate               (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_resolution_changed                  (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_bgcolor_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_bg_color_selector_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_portable_pixmap1_activate           (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_postscript1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_encapsulated_postscript1_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_portable_document_format1_activate  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_latex1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_none1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_barycenter1_activate                (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_bsp_tree1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_portrait1_activate                  (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_landscape1_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_line_offset_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_best_root_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_no_text_toggled                     (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_no_pixmap_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_no_ps3_toggled                      (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_browse_clicked                      (GtkButton       *button,
                                        gpointer         user_data);

void
on_save_as_button_clicked              (GtkButton       *button,
                                        gpointer         user_data);

void
on_lines_closer_changed                (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_about1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_front1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_delete_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_entry2_changed                      (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_gnuplot1_activate                   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_maxlevel_changed                    (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_reactivity_changed                  (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_spinbuttonx_changed                 (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_spinbuttony_changed                 (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_spinbuttonz_changed                 (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_spinbuttonpos_changed               (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_spinbuttonlevel_changed             (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_reversed_toggled                    (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_solid_reversed_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

gboolean
main_quit                              (GtkWidget       *widget,
                                        GdkEvent        *event,
                                        gpointer         user_data);

void
on_spinbuttoniso_value_changed         (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_back1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_top1_activate                       (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_bottom1_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_left1_activate                      (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_right1_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_ellipse_scale_changed               (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_symbol_size_value_changed           (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_substract1_activate                 (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_back_clicked                        (GtkButton       *button,
                                        gpointer         user_data);

void
on_forward_clicked                     (GtkButton       *button,
                                        gpointer         user_data);

void
on_gfsview_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_gerris_activate                     (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_scripting_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_entryiso_activate                   (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_scalar_changed                      (GtkEditable     *editable,
                                        gpointer         user_data);

void
mark_edited                            (GtkEditable     *editable,
                                        gpointer         user_data);

void
unmark_edited                          (GObject          *object,
                                        gpointer         user_data);

void
on_finest_toggled                      (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_delete_streamline_clicked           (GtkButton       *button,
                                        gpointer         user_data);

void
on_new_streamline_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_show_cells_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_dmin_changed                        (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_evenly_spaced_clicked               (GtkButton       *button,
                                        gpointer         user_data);

void
on_stop_evenly_clicked                 (GtkButton       *button,
                                        gpointer         user_data);

void
on_radius_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

gboolean 
gl2ps_update_ppm_size                  (GtkWidget * widget, 
					GdkEventExpose * event, 
					GtkWidget * gl2ps);

void
on_screen_ppm_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_width_value_changed                 (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_height_value_changed                (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_keep_ratio_toggled                  (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
mark_edited                            (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_information_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_locate_x_value_changed              (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_locate_y_value_changed              (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_locate_z_value_changed              (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_locate_x_entry_changed              (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_locate_y_entry_changed              (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_locate_z_entry_changed              (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_x_scaling_changed                   (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_y_scaling_changed                   (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_z_scaling_changed                   (GtkEditable     *editable,
                                        gpointer         user_data);

void
on_vof_reversed_toggled                (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_line_width_changed                  (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_scalable_vector_graphics_activate   (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_vof_visible_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_wavefront_obj1_activate             (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_kml_obj_activate                    (GtkMenuItem     *menuitem,
                                        gpointer         user_data);

void
on_linear_reversed_toggled             (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_location_label_toggled              (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_font_size_changed                   (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_vector_font_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_raster_font_toggled                 (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_labelx_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_labely_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_labelz_changed                      (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);

void
on_labelentry_activate                 (GtkEntry        *entry,
                                        gpointer         user_data);

void
on_label_symbol_check_toggled          (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_showscale_toggled                   (GtkToggleButton *togglebutton,
                                        gpointer         user_data);

void
on_linewidth_value_changed             (GtkSpinButton   *spinbutton,
                                        gpointer         user_data);
