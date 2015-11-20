#include "consolidate.h"
#include "diskreadmda.h"
#include "mdaio.h"

void get_load_channels(int K,int *load_channels,DiskReadMda &cluster);

bool consolidate(const char *cluster_path,const char *templates_path,const char *cluster_out_path,const char *templates_out_path,const char *load_channels_out_path) {
    DiskReadMda cluster; cluster.setPath(cluster_path);
    DiskReadMda templates; templates.setPath(templates_path);

    int M=templates.N1();
    int T=templates.N2();
    int K=templates.N3();

    int load_channels[K+1];
    get_load_channels(K,load_channels,cluster);

    bool to_use[K+1];
    int num_to_use=0;
    for (int k=1; k<=K; k++) {
        float energies[M+1];
        for (int m=1; m<=M; m++) {
            float sumsqr=0;
            for (int t=0; t<T; t++) {
                float val=templates.value(m-1,t,k-1);
                sumsqr+=val*val;
            }
            energies[m]=sumsqr;
        }
        int ch=load_channels[k];
        bool okay=true;
        for (int m=1; m<=M; m++) {
            if (m!=ch) {
                if (energies[ch]<=0.9*energies[m]) {
                    okay=false;
                }
            }
        }
        to_use[k]=okay;
        if (okay) num_to_use++;
    }

    int K2=num_to_use;
    int num_events2=0;
    for (int i=0; i<cluster.N2(); i++) {
        int k=(int)cluster.value(2,i);
        if (to_use[k]) num_events2++;
    }

    int kk=1;
    int cluster_mapping[K+1];
    for (int i=1; i<=K; i++) {
        if (to_use[i]) {
            cluster_mapping[i]=kk;
            kk++;
        }
        else {
            cluster_mapping[i]=0;
        }
    }

    MDAIO_HEADER H_cluster;
    H_cluster.data_type=MDAIO_TYPE_FLOAT32;
    H_cluster.num_bytes_per_entry=4;
    H_cluster.num_dims=2;
    H_cluster.dims[0]=3;
    H_cluster.dims[1]=num_events2;
    {
        FILE *outf=fopen(cluster_out_path,"wb");
        if (!outf) {
            printf("Unable to open file for writing: %s\n",cluster_out_path);
            return false;
        }
        mda_write_header(&H_cluster,outf);
        for (int i=0; i<cluster.N2(); i++) {
            int k0=(int)cluster.value(2,i);
            if (to_use[k0]) {
                int k=cluster_mapping[k0];
                float ch_tt_k[3];
                ch_tt_k[0]=cluster.value(0,i);
                ch_tt_k[1]=cluster.value(1,i);
                ch_tt_k[2]=k;
                mda_write_float32(ch_tt_k,&H_cluster,3,outf);
            }
        }
        fclose(outf);
    }


    MDAIO_HEADER H_templates;
    H_templates.data_type=MDAIO_TYPE_FLOAT32;
    H_templates.num_bytes_per_entry=4;
    H_templates.num_dims=3;
    H_templates.dims[0]=M;
    H_templates.dims[1]=T;
    H_templates.dims[2]=K2;
    {
        FILE *outf=fopen(templates_out_path,"wb");
        if (!outf) {
            printf("Unable to open file for writing: %s\n",templates_out_path);
            return false;
        }
        mda_write_header(&H_templates,outf);
        for (int i=1; i<=K; i++) {
            if (to_use[i]) {
                float *buf=(float *)malloc(sizeof(float)*M*T);
                int aa=0;
                for (int t=0; t<T; t++) {
                    for (int m=0; m<M; m++) {
                        buf[aa]=templates.value(m,t,i-1);
                        aa++;
                    }
                }
                mda_write_float32(buf,&H_templates,M*T,outf);
                free(buf);
            }
        }
        fclose(outf);
    }

    MDAIO_HEADER H_load_channels;
    H_load_channels.data_type=MDAIO_TYPE_FLOAT32;
    H_load_channels.num_bytes_per_entry=4;
    H_load_channels.num_dims=2;
    H_load_channels.dims[0]=1;
    H_load_channels.dims[1]=K2;
    {
        FILE *outf=fopen(load_channels_out_path,"wb");
        if (!outf) {
            printf("Unable to open file for writing: %s\n",load_channels_out_path);
            return false;
        }
        mda_write_header(&H_load_channels,outf);
        for (int i=1; i<=K; i++) {
            if (to_use[i]) {
                float val=load_channels[i];
                mda_write_float32(&val,&H_load_channels,1,outf);
            }
        }
        fclose(outf);
    }

    return true;
}

void get_load_channels(int K,int *load_channels,DiskReadMda &cluster) {
    for (int k=0; k<=K; k++) load_channels[k]=0;
    for (int i=0; i<cluster.N2(); i++) {
        int ch=(int)cluster.value(0,i);
        int k=(int)cluster.value(2,i);
        if ((1<=k)&&(k<=K)) load_channels[k]=ch;
    }
}
