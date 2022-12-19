from django.views import generic
from accounts.models import Profile
from django.contrib.auth.models import User
from blog.models import Comments, Posts, Category
from hitcount.models import Hit
from django.contrib import messages
from django.shortcuts import redirect
from django.conf import settings
from blog.models import Comments
from django_plotly_dash import DjangoDash, apps
from dash import Dash, html, dcc, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
from collections import OrderedDict


class DashboardView(generic.ListView):
    template_name = "dashboard/pages/dashboard.html"
    model = Profile

    def get(self, request, *args, **kwargs):
        if not request.user.is_staff:
            messages.error(request, "Yetkili girişi yapınız!")
            return redirect('%s?next=/blog/' % (settings.LOGIN_URL))

        user_agent = []
        users = []
        created = []

        for hit in Hit.objects.all():
            user_agent.append(hit.user_agent)
            created.append(hit.created.date())
            users.append(str(hit.user))

        hits_table = DjangoDash('Hits')

        data = OrderedDict(
            [
                ("Tarih", created),
                ("Kullanıcılar", users),
                ("İP", user_agent),
            ]
        )

        hits_data = pd.DataFrame(
            OrderedDict(
                [(name, col_data) for (name, col_data) in data.items()]
            )
        )

        hits_table.layout = html.Div([
            dbc.Label('Kullanıcı - IP'),
            dash_table.DataTable(
                data=hits_data.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in hits_data.columns],
                page_size=10,
                style_table={'height': '250px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
            )
        ])

        user_graph = DjangoDash('DateUser')

        users_table = DjangoDash('UsersTable')

        date_joined = []
        user_model = []
        first_name = []
        last_login = []
        last_login_time = []
        last_name = []

        for user in User.objects.all():
            user_model.append(user.username)
            first_name.append(user.first_name)
            last_name.append(user.last_name)
            last_login_time.append(user.last_login.time())
            last_login.append(user.last_login.date())
            date_joined.append(user.date_joined.date())

        users_data = OrderedDict(
            [
                ("Katılma Tarihi", date_joined),
                ("Kullanıcılar", user_model),
                ("Adı", first_name),
                ("Soyadı", last_name),
                ("Son görülme", last_login),
                ("Son görülme saat", last_login_time),
            ]
        )

        users_data_frame = pd.DataFrame(
            OrderedDict([(name, col_data) for (name, col_data) in users_data.items()])
        )

        users_table.layout = html.Div([
            dbc.Label('Kullanıcı'),
            dash_table.DataTable(
                data=users_data_frame.to_dict('records'),
                columns=[{'id': c, 'name': c} for c in users_data_frame.columns],
                page_size=10,
                style_table={'height': '250px', 'overflowY': 'auto'},
                style_cell={'textAlign': 'center'},
            )
        ])

        fig = px.line(data_frame=users_data_frame, x="Katılma Tarihi")

        user_graph.layout = html.Div([
            dbc.Label('Kullanıcılar'),
            dcc.Graph(figure=fig)
        ])

        return super().get(request, *args, **kwargs)

    def get_context_data(self, *, object_list=None, **kwargs):
        context = super().get_context_data(**kwargs)
        context['comments'] = Comments.objects.filter(post__author=self.request.user).order_by('-id')[:10]
        return context
