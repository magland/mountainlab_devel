#include "mvclipsview.h"

class MVClipsViewPrivate {
public:
    MVClipsView *q;
    Mda m_clips;
    QList<double> m_times;
};

MVClipsView::MVClipsView()
{
    d=new MVClipsViewPrivate;
    d->q=this;
    this->initialize();
    connect(this,SIGNAL(currentXChanged()),this,SIGNAL(currentClipTimepointChanged()));
}

MVClipsView::~MVClipsView()
{
    delete d;
}

void MVClipsView::setClips(const Mda &clips)
{
    d->m_clips=clips;
    DiskArrayModel *DAM=new DiskArrayModel;
    DAM->setFromMda(d->m_clips);
    this->setData(DAM,true);
}

void MVClipsView::setTimes(const QList<double> &times)
{
    d->m_times=times;
}

int MVClipsView::currentClipIndex()
{
    int T=d->m_clips.N2();
    double tp=this->currentX();
    int clip_index=(int)(tp/T);
    return clip_index;
}

double MVClipsView::currentClipTimepoint()
{
    int T=d->m_clips.N2();
    int Tmid=(int)((T+1)/2)-1;
    double tp=this->currentX();
    int clip_index=(int)(tp/T);
    double offset=tp-(clip_index*T)-Tmid;
    qDebug() << tp << clip_index << offset;
    return d->m_times.value(clip_index)+offset;
}

