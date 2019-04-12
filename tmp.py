custom_hover = []
  for _, i in neut_df.iterrows():
      dat = 'Gene: '+'{}'.format(i['Gene'])+'<br>'+'Feature: '+'{}'.format(i['Feature'])+'<br>'+'Value: '+'{:.2f}'.format(i['Value'])+'<br>'+'R: '+'{:.2f}'.format(i['R'])
      custom_hover.append(dat)

  return {
      'data': [(
          go.Heatmap(
              x=neut_df['Gene'],
              y=neut_df['Feature'],
              z=neut_df['Value'],
              name='neut-heatmap',
              colorscale='RdBu',
              text=custom_hover,
              hoverinfo='text')
              )],
      'layout':go.Layout(
          title=go.layout.Title(
              text=('<b>Target label: '+neut_df['Type'].iloc[0]+'</b>'),
              xanchor='left',
              yanchor='bottom',
              x=0,
              font=dict(
                  family='Arial',
                  size=16,
                  color='black'
              )
          ),
          autosize=False,
          yaxis=dict(
              automargin=True,
              tickfont=dict(
                  family='Arial, sans-serif',
                  size=14,
                  color='black'
              )
          )
      )
  }
