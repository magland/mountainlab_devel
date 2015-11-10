#ifndef FIRETRACKWIDGET_H
#define FIRETRACKWIDGET_H

#include <QMainWindow>
#include <QWheelEvent>
#include "diskreadmda.h"
#include "mda.h"

class FireTrackWidgetPrivate;

class FireTrackWidget : public QMainWindow
{
	Q_OBJECT
public:
	friend class FireTrackWidgetPrivate;
	explicit FireTrackWidget(QWidget *parent = 0);

	void setElectrodeLocations(const Mda &L);
	void setWaveforms(const Mda &X);

    Q_INVOKABLE void animate();
	Q_INVOKABLE void setWaveforms(DiskReadMda *X);
	Q_INVOKABLE void setElectrodeLocations(DiskReadMda *X);

protected:
	void resizeEvent(QResizeEvent *evt);

private slots:
	void slot_current_waveform_changed();
	void slot_selected_electrodes_changed();
	void slot_timepoint_changed();
	void slot_plot_timepoint_changed();
	void slot_electrode_left_clicked(int ind);
	void slot_plot_options_changed();
	void slot_plot_vertical_scaling(float val);
    void slot_loop();

private:
	FireTrackWidgetPrivate *d;

	~FireTrackWidget();

signals:

public slots:
};

#endif // FIRETRACKWIDGET_H
